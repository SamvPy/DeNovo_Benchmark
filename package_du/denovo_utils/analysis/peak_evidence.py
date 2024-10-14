"""Utilities to extract peak evidence for a PSM."""

import warnings

import numpy as np
from ms2pip.result import ProcessingResult
from pyteomics.mass.mass import Composition

warnings.filterwarnings("ignore")

LOG_NULL = np.float32(np.round(np.log2(0.001), 6))

UNIMOD_MOD_COMP = {
    "[UNIMOD:1]": Composition({"H": 2, "C": 2, "O": 1}),
    "[UNIMOD:4]": Composition({"H": 3, "C": 2, "N": 1, "O": 1}),
    "[UNIMOD:5]": Composition({"H": 1, "C": 1, "N": 1, "O": 1}),
    "[UNIMOD:7]": Composition({"H": -1, "N": -1, "O": 1}),
    "[UNIMOD:21]": Composition({"H": 1, "O": 3, "P": 1}),
    "[UNIMOD:23]": Composition({"H": -2, "O": -1}),
    "[UNIMOD:35]": Composition({"O": 1}),
    "[UNIMOD:385]": Composition({"H": -3, "N": -1}),
    "[+25.980265]": Composition({"H": -2, "C": 1, "O": 1}),
}
UNIMOD_LIST = [
    "[UNIMOD:1]",
    "[UNIMOD:4]",
    "[UNIMOD:5]",
    "[UNIMOD:7]",
    "[UNIMOD:21]",
    "[UNIMOD:23]",
    "[UNIMOD:35]",
    "[UNIMOD:385]",
]


class Evidence:
    """Evidence class providing tools to see supporting ions of a PSM."""

    def __init__(self, processing_result):
        """
        Initialize peak evidence object for a peptide sequence.

        Can be used to easily compare peptide sequences while allowing
        isobaric segments.

        Parameter
        ---------
        spectrum: dict or ProcessingResult
            A dictionary with keys 'observed_intensity',
            'theoretical_mz' and 'peptide
        """
        self.peptide = processing_result.psm.peptidoform.proforma
        self.processing_result = processing_result
        self.evidence = self.get_evidence(processing_result=processing_result)
        # self.peptide_masked = self.mask_peptide()

    def get_evidence(
        self, processing_result: ProcessingResult, ion_types=["b", "b2", "y", "y2"]
    ):
        """
        Get peak evidence for a PSM.

        The matching is done with annotate_spectra from ms2pip package.
        To have evidence of an amino acid, I require a fragment ion series.

        (Could include an implementation where the matching
        is performed with other evidence peaks as well...)
        """
        psm = processing_result.psm
        _ = processing_result.theoretical_mz
        intensities = processing_result.observed_intensity

        # Ensure the floats are of the same type
        intensities = {
            k: np.array([np.round(i, 6) for i in v], dtype=np.float32)
            for k, v in intensities.items()
        }
        observed_bool = {k: v != LOG_NULL for k, v in intensities.items()}

        evidence = []
        residue_list = [aa for aa in psm.peptidoform]

        for i, _ in enumerate(psm.peptidoform):
            evidence_types = {}
            evidence_present = False

            if i == len(psm.peptidoform) - 1:
                break

            for ion_type in ion_types:
                if peak_is_observed(
                    observed_list=observed_bool[ion_type],
                    ion_type=ion_type,
                    ion_number=i,
                ):
                    evidence_present = True
                    evidence_types[ion_type] = True
                else:
                    evidence_types[ion_type] = False

            evidence.append(
                {
                    "evidence": evidence_present,
                    "evidence_types": evidence_types,
                    "cleavage_site_residues": [residue_list[i], residue_list[i + 1]],
                    "cleavage_site_indices": [i, i + 1],
                }
            )

        return evidence

    def mask_peptide(self):
        """To be implemented."""
        pass
        # self.evidence

    # def __repr__(self):
    #     return self.evidence

    def __eq__(self):
        """
        Allow isobaric alternatives in a portion of the sequence without peak evidence.

        To be implemented.
        """
        pass


def peak_is_observed(observed_list, ion_type, ion_number):
    """Evaluate if theoretical peak is observed."""
    # Y peak evidence for amino acid 1 is the yn peak!
    if "y" in ion_type:
        observed_list = observed_list[::-1]
    return observed_list[ion_number]
