"""Class for feature generation."""

from typing import List

import pyopenms as oms
from ms2rescore.feature_generators.base import FeatureGeneratorBase
from psm_utils import PSM, PSMList
from pyteomics.mgf import IndexedMGF

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
