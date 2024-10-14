from psm_utils import PSM
import logging
import re
from typing import Literal, Optional

from tqdm import tqdm

import pyopenms as oms
from psm_utils import Peptidoform
from psm_utils.peptidoform import PeptidoformException
from pyteomics.proforma import ProFormaError
from spectrum_utils import proforma
from spectrum_utils.fragment_annotation import get_theoretical_fragments

from ms2rescore.feature_generators.base import FeatureGeneratorBase
from typing import List
from psm_utils import PSM, PSMList

def proforma_to_oms(peptide: str) -> tuple[oms.AASequence, Optional[int]]:
    """
    Parse a peptide sequence in proforma format to pyOpenMS compatible format.

    Parameter
    ---------
    peptide: str
        Peptide string in proforma format

    Returns
    -------
    AASequence (pyOpenMS):
        A peptide sequence in pyOpenMS format
    int:
        charge of the peptide
    """
    # Check if peptide is in proforma format
    try:
        _ = Peptidoform(peptide)
    except ProFormaError:
        raise PeptidoformException(f"({peptide}) not in proforma format.")

    # Error with UNIMOD:385
    peptide = peptide.replace("UNIMOD:385", "-17.027")
    peptide = peptide.replace("UNIMOD:5", "43.005814")
    peptide = peptide.replace("Formula:H-2C1O1", "+25.980265")

    # Reformat unimod modifications
    pattern_unimod = r"\[UNIMOD:(\d+)\]"

    def replace_unimod(match):
        return f"(UniMod:{match.group(1)})"

    peptide_oms_str = re.sub(
        pattern=pattern_unimod, repl=replace_unimod, string=peptide
    )

    # Parse N-terminal modifications
    if ")-" in peptide_oms_str:
        peptide_oms_list = peptide_oms_str.split(")-")
        nterm_modification, peptide_oms_str = peptide_oms_list[-2], peptide_oms_list[-1]
        nterm_modification += ")"
        peptide_oms_str = "." + nterm_modification + peptide_oms_str + "."
    elif "]-" in peptide_oms_str:
        peptide_oms_list = peptide_oms_str.split("]-")
        nterm_modification, peptide_oms_str = peptide_oms_list[-2], peptide_oms_list[-1]
        nterm_modification += "]"
        peptide_oms_str = "." + nterm_modification + peptide_oms_str + "."

    # Split the charge from the peptide string
    if "/" in peptide_oms_str:
        peptide_oms_str, charge = peptide_oms_str.split("/")
    else:
        charge = None

    peptide_oms = oms.AASequence.fromString(peptide_oms_str)

    return peptide_oms, charge


def proforma_to_theoretical_spectrum(
    peptide: str, engine: Literal["spectrum-utils", "pyopenms"] = "spectrum-utils"
) -> oms.MSSpectrum:
    """
    Create a theoretical spectrum from a peptide sequence.

    Parameter
    ---------
    peptide: str
        Peptide sequence in proforma format
    engine: str
        The engine to use to create theoretical spectrum.
        Can only be 'pyopenms' or 'spectrum-utils' (default)

    Return
    ------
    MSSpectrum
        Spectrum object in pyOpenMS format
    """
    if engine == "spectrum-utils":
        proforma_str = proforma.parse(peptide)[0]
        theoretical_fragments = get_theoretical_fragments(proteoform=proforma_str)
        mz_array = [peak[1] for peak in theoretical_fragments]
        intensity_array = [1.0] * len(mz_array)

        spectrum = oms.MSSpectrum()
        spectrum.set_peaks([mz_array, intensity_array])
        return spectrum

    elif engine == "pyopenms":
        # Reformat peptide sequence in pyOpenMS format
        peptide_oms, charge = proforma_to_oms(peptide=peptide)

        # Initialize the required objects to create the spectrum
        spectrum = oms.MSSpectrum()
        tsg = oms.TheoreticalSpectrumGenerator()
        p = oms.Param()

        p.setValue("add_b_ions", "true")
        p.setValue("add_metainfo", "true")
        tsg.setParameters(param=p)

        # Create the theoretical spectrum
        tsg.getSpectrum(spec=spectrum, peptide=peptide_oms, min_charge=1, max_charge=2)
        return spectrum

    else:
        raise Exception("Engine not suported for theoretical spectrum generation.")


def calculate_hyperscore(
    psm: PSM,
    spectrum: dict,
    fragment_tol_mass=20,
    fragment_tol_mode="ppm"
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
        peptide=psm.peptidoform.proforma, engine="pyopenms"
    )
    observed_spectrum_oms = oms.MSSpectrum()
    observed_spectrum_oms.set_peaks(
        [spectrum["m/z array"], spectrum["intensity array"]]
    )
    hyperscore = oms.HyperScore()
    result = hyperscore.compute(
        fragment_mass_tolerance=fragment_tol_mass,
        fragment_mass_tolerance_unit_ppm=fragment_mass_tolerance_unit_ppm,
        exp_spectrum=observed_spectrum_oms,
        theo_spectrum=theoretical_spectrum,
    )
    return result

class HyperscoreGenerator(FeatureGeneratorBase):
    """MS2Rescore type feature generator for ppm-errors"""
    def __init__(self, *args, **kwargs):
        """Initialize feature generator class."""
        self.config = kwargs
        super().__init__(*args, **kwargs)

    @property
    def feature_names(self) -> List[str]:
        """Names of features added to rescoring_features dict."""
        return [
            "hyperscore"
        ]

    @property
    def name(self) -> str:
        return "hyperscore"
    @property
    def input_type(self) -> str:
        return "psm_list"
    
    def add_features(self, psm_list: PSMList, spectra: List) -> None:
        """Compute and add rescoring features to psmlist."""

        for psm, spectrum in tqdm(zip(psm_list, spectra)):
            psm.rescoring_features.update({"hyperscore": calculate_hyperscore(
                    psm=psm,
                    spectrum=spectrum,
                    fragment_tol_mass=self.config["fragment_tol_mass"],
                    fragment_tol_mode=self.config["fragment_tol_mode"] 
                )
            })
