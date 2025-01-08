from pyteomics import mgf
import polars as pl
import os
from glob import glob
import numpy as np
from rustyms import RawSpectrum, LinearPeptide, FragmentationModel

from ..utils.utils import (
    get_annotated_spectrum,
    annot_peaks_to_fragments,
    fragments_to_polars,
    parse_to_iontype_dict,
    ion_dict_to_matrix,
    matrix_to_ion_dict,
    calculate_ppm,
    calculate_ppm_matrix
)

class SpectrumVector():
    def __init__(
            self,
            ion_types,
            neutral_losses,
    ):
        self.annotated = False
        self.ion_types = ion_types
        self.neutral_losses = neutral_losses
        self.tic = 0.0
        self.precursor_mz = np.nan

    def parse(self, psm, spectrum, spectrum_format='mgf'):
        annot_spec, theo_frags = get_annotated_spectrum(psm, spectrum)
        annot_spec = annot_spec.spectrum
        self.by = get_by_fragments(annot_spec)
        self.tic = np.sum(spectrum["intensity array"])
        self.peptidoform = psm.peptidoform
        self.spectrum_id = psm.spectrum_id

        if spectrum_format=='mgf':
            self.precursor_mz = spectrum["params"]["pepmass"][0]
        elif spectrum_format=='mzML':
            self.precursor_mz = spectrum['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]['selected ion m/z']
        else:
            raise Exception(f"Selected spectrum format ({spectrum_format}) is not supported in SpectrumVector parsing.")

        self.precursor_ppm = calculate_ppm(
            # Observed mz (extracted from raw file)
            m1=self.precursor_mz,
            # Theoretical mz, which will be modified with isotopes
            m2=psm.peptidoform.theoretical_mz,
            charge=psm.peptidoform.precursor_charge,
            isotopes=[-1,0,1,2],
            max_value=50
        )

        # Convert the list of fragments to a polars object
        annot_frags, mz_array, intensity_array = annot_peaks_to_fragments(annot_spec)

        spec_polars_theo = fragments_to_polars(
            fragment_list=theo_frags,
            ion_types=self.ion_types,
            neutral_losses=self.neutral_losses
        )
        spec_polars_annot = fragments_to_polars(
            fragment_list=annot_frags,
            ion_types=self.ion_types,
            neutral_losses=self.neutral_losses,
            mz_array=mz_array,
            intensity_array=intensity_array
        )
        
        # Parse the polars spectrum to a dictionary
        result_dict_theo_mz = parse_to_iontype_dict(
            pl_df=spec_polars_theo,
            len_pep=len(psm.peptidoform),
            ion_types=self.ion_types,
            value="mz"
        )
        result_dict_annot_mz = parse_to_iontype_dict(
            pl_df=spec_polars_annot,
            len_pep=len(psm.peptidoform),
            ion_types=self.ion_types,
            value="mz"
        )
        result_dict_annot_intensity = parse_to_iontype_dict(
            pl_df=spec_polars_annot,
            len_pep=len(psm.peptidoform),
            ion_types=self.ion_types,
            value="intensity"
        )

        self.theoretical_mz = result_dict_theo_mz
        self.experimental_mz = self._complete_ion_dict(result_dict_annot_mz)
        self.experimental_intensity = self._complete_ion_dict(result_dict_annot_intensity)
        self.annotated=True

    def _complete_ion_dict(self, ion_dict):
        for nl in self.neutral_losses:
            if nl not in ion_dict.keys():
                ion_dict[nl] = {i: np.zeros(self.n) for i in self.ion_types}
            else:
                for ion_type in self.ion_types:
                    if ion_type not in ion_dict[nl].keys():
                        ion_dict[nl][ion_type] = np.zeros(self.n)
        return ion_dict

    def load(
            self,
            theoretical_mz,
            experimental_mz,
            experimental_intensity,
            tic,
            peptidoform,
            spectrum_id,
            precursor_mz,
            precursor_ppm
        ):
        """
        Load a parsed spectrum.
        """
        self.annotated = True
        self.theoretical_mz = theoretical_mz
        self.experimental_mz = experimental_mz
        self.experimental_intensity = experimental_intensity
        self.tic = tic
        self.peptidoform = peptidoform
        self.spectrum_id = spectrum_id
        self.precursor_mz = precursor_mz
        self.precursor_ppm = precursor_ppm

    def delete_shared_peaks(self, method):
        self.experimental_mz, self.experimental_intensity = method(self.experimental_mz, self.experimental_intensity)

    @property
    def ppm_diff(self):
        matrix_theoretical = ion_dict_to_matrix(
            self.theoretical_mz[""],
            ion_types=self.ion_types,
            n=len(self.theoretical_mz[""]["y1"])
        )
        matrix_theoretical[matrix_theoretical==0] = 1
        matrix_experimental = ion_dict_to_matrix(
            self.experimental_mz[""],
            ion_types=self.ion_types,
            n=len(self.theoretical_mz[""]["y1"])
        )

        ppm_diff = np.where(
            matrix_experimental!=0, 
            calculate_ppm_matrix(
                matrix_theoretical,
                matrix_experimental,
            ),
            np.nan
        )
        return matrix_to_ion_dict(ppm_diff, ion_types=self.ion_types)

    @property
    def n(self):
        return len(self.peptidoform)-1
    

def infer_fragment_identity(frag, allow_ion_types=['b', 'y']):
    ion = frag.ion

    is_allowed = False
    for allowed_ion_type in allow_ion_types:
        if allowed_ion_type in ion:
            is_allowed=True
            break
    
    if not is_allowed:
        return False
    # Radicals
    if "·" in ion:
        return False
    if frag.neutral_loss is not None:
        return False
    if frag.charge > 2:
        return False
    
    return ion[0]

def get_by_fragments(annotated_spectrum):
    b_intensities = []
    y_intensities = []
    for peak in annotated_spectrum:
        
        for fragment in peak.annotation:

            ion_type = infer_fragment_identity(fragment)
            
            if ion_type == 'b':
                b_intensities.append(peak.intensity)
            if ion_type == 'y':
                y_intensities.append(peak.intensity)
    return b_intensities, y_intensities