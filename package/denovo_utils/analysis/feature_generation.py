"""Class for feature generation."""

from typing import List

import pyopenms as oms
from ms2rescore.feature_generators.base import FeatureGeneratorBase
from psm_utils import PSM, PSMList
from pyteomics.mgf import IndexedMGF
import numpy as np

from ..parsers import proforma_to_theoretical_spectrum


class PeakFeatureGenerator(FeatureGeneratorBase):
    """MS2Rescore type feature generator for peak-type features."""

    def __init__(self, *args, **kwargs):
        """Initialize feature generator class."""
        super().__init__(*args, **kwargs)
        NotImplementedError()

    @property
    def feature_names(self) -> List[str]:
        """Names of features added to rescoring_features dict."""
        return ["hyperscore", "peak_tic", "peak_count"]

    def add_features(self, psm_list: PSMList) -> None:
        """Compute and add rescoring features to psmlist."""
        NotImplementedError()


class ExplainedIntensityFeatures():
    def __init__(self):
        self.features = {}

    def get_sum_intensity_iontype(self, spectrum, ion_type):
        spectrum.isotope_matrix[ion_type]

    def add_features(self, spectrum):
        # Define a noise cutoff
        annotation_mask = get_annotated_mask(spectrum)
        noise_i, noise_mz = spectrum.spectrum.intensity[~annotation_mask], spectrum.spectrum.mz[~annotation_mask]
        intensity_explained = np.sum(spectrum.spectrum.intensity[annotation_mask])
        
        noise_cutoff = set_noise_cutoff(noise_i=noise_i, lod_type="lod_10")
        above_noise = spectrum.spectrum.intensity > noise_cutoff

        intensity_explained_above_noise = np.sum(spectrum.spectrum.intensity[np.logical_and(
            annotation_mask, above_noise)
        ])

        self.features["pct_peaks_explained"] = np.sum(annotation_mask) / len(annotation_mask)
        self.features["pct_intensity_explained"] = intensity_explained / spectrum.tic
        
        if np.sum(above_noise)==0:
            self.features["pct_peaks_explained_above_noise"]=0
            self.features["pct_intensity_explained_above_noise"]=0

        else:
            self.features["pct_peaks_explained_above_noise"] = np.sum(np.logical_and(
                above_noise,
                annotation_mask
            )) / np.sum(above_noise)
            self.features["pct_intensity_explained_above_noise"] = intensity_explained_above_noise / np.sum(spectrum.spectrum.intensity[above_noise])

        self.features["n_big_unannotated_peaks"] = np.sum(np.logical_and(
            above_noise, ~annotation_mask
        ))
        
        for ion_type in "y1 y2 b1 b2 x1 x2 a1 a2 c1 c2 z1 z2".split():
            self.features[f"{ion_type}_pct_explained"] = np.sum(spectrum.isotope_matrix[ion_type]) / intensity_explained
        self.features["precursor_pct_explained"] = np.sum(spectrum.isotope_matrix["p"]) / intensity_explained

        by_intensity = sum([np.sum(spectrum.isotope_matrix[i]) for i in "y1 y2 b1 b2".split()])
        self.features["pct_precursor_fragmented"] = by_intensity / (by_intensity + np.sum(spectrum.isotope_matrix["p"]))
        self.features["noise_cutoff"] = noise_cutoff


def set_noise_cutoff(noise_i, lod_type="lod_5"):
    """Select noise cutoff
    
    Options:
    - lod_5
    - lod_10
    - std_cutoff"""
    
    if lod_type == "std_cutoff":
        log_i = np.log2(noise_i)
        log_cutoff = np.median(log_i)+np.std(log_i)
        int_cutoff = 2**log_cutoff
        return int_cutoff

    elif lod_type == "lod_5":
        lod_5 = min(noise_i)*5
        return lod_5

    elif lod_type == "lod_10":
        lod_10 = min(noise_i)*10
        return lod_10
    
    else:
        NotImplementedError()


def calculate_hyperscore(
    psm: PSM,
    mgf_file: IndexedMGF,
    fragment_tol_mass=10,
    fragment_tol_mode="ppm",
    engine="spectrum-utils",
):
    """
    Calculate the hyperscore as defined in the X!Tandem search engine.

    It is a metric of how good two spectra match with each other (matching peaks).

    Parameters
    ----------
    psm: psm_utils.PSM
        The PSM used to extract 'spectrum_id' (for MGF spectrum extraction)
        and 'Peptidoform' (the peptide sequence)
    mgf_file: pyteomics.mgf.IndexedMGF
        An mgf-file (or mzml) used to extract the spectrum
        (based on 'spectrum_id' in PSM and 'id' in the file)
    fragment_tol_mass: int
        The allowed tolerance to match peaks
    fragment_tol_mode: str
        'ppm' for parts-per-million mode. 'Da' for fragment_tol_mass in Dalton.
    engine: str (default: spectrum-utils)
        Engine to use for creating theoretical spectrum. Defaults to spectrum-utils.
        Other options are 'pyopenms'.

    Return
    ------
    int
        The hyperscore
    """
    if fragment_tol_mode == "ppm":
        fragment_mass_tolerance_unit_ppm = True
    elif fragment_tol_mode == "Da":
        fragment_mass_tolerance_unit_ppm = False
    else:
        raise Exception(
            "fragment_tol_mode can only take 'Da' or 'ppm'. {} was provided.".format(
                fragment_tol_mode
            )
        )

    theoretical_spectrum = proforma_to_theoretical_spectrum(
        peptide=psm.peptidoform.proforma, engine=engine
    )
    observed_spectrum = mgf_file.get_by_id(psm.spectrum_id)
    observed_spectrum_oms = oms.MSSpectrum()
    observed_spectrum_oms.set_peaks(
        [observed_spectrum["m/z array"], observed_spectrum["intensity array"]]
    )

    hyperscore = oms.HyperScore()
    result = hyperscore.compute(
        fragment_mass_tolerance=fragment_tol_mass,
        fragment_mass_tolerance_unit_ppm=fragment_mass_tolerance_unit_ppm,
        exp_spectrum=observed_spectrum_oms,
        theo_spectrum=theoretical_spectrum,
    )
    return result
