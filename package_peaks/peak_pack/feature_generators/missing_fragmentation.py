from ms2rescore.feature_generators.base import FeatureGeneratorBase
from typing import List
from psm_utils import PSM, PSMList

import numpy as np
from tqdm import tqdm

from ..annotation.evidence import PeptideEvidence


class MissingFragmentationFeatures(FeatureGeneratorBase):
    """MS2Rescore type feature generator for peak-type features."""

    def __init__(self, *args, **kwargs):
        """Initialize feature generator class."""
        super().__init__(*args, **kwargs)

    @property
    def feature_names(self) -> List[str]:
        """Names of features added to rescoring_features dict."""
        return [
            "missing_frag_count",
            "missing_frag_pct",
            "missing_frag_sites",
            "missing_frag_longest_sequence",
            "missing_frag_from_n",
            "missing_frag_from_c",
        ]

    @property
    def name(self) -> str:
        return "missing_frag"
    @property
    def input_type(self) -> str:
        return "peptide_evidence"

    def add_features(self, psm_list: PSMList, peptide_evidence_list: List[PeptideEvidence]) -> None:
        for psm, pe in tqdm(zip(psm_list, peptide_evidence_list)):

            psm.rescoring_features.update(
                self.calculate_features(pe.evidence)
            )

    def calculate_features(self, evidence):
        # Ensure it's a boolean array
        arr = np.asarray(evidence, dtype=bool)
        
        # Number of consecutive False entries starting from the start
        consecutive_false_start = np.argmax(arr) if not arr.all() else 0

        # Number of consecutive False entries starting from the end
        consecutive_false_end = np.argmax(arr[::-1]) if not arr.all() else 0

        # Number of False entries in the array
        num_false_entries = np.size(arr) - np.sum(arr)

        # Percentage of False entries in the array
        percentage_false = num_false_entries / len(arr) *100

        # Finding the longest sequence of consecutive False values
        # Identify the regions where values change (False <-> True)
        # Create an array to mark transitions
        if not arr.all():
            # Extend the array with True at both ends to handle edge cases
            padded_arr = np.r_[True, arr, True]
            diff_arr = np.diff(padded_arr.astype(int)) # Transitions
            false_starts = np.flatnonzero(diff_arr == -1)
            false_ends = np.flatnonzero(diff_arr == 1)
            longest_false_seq = (false_ends - false_starts).max() if len(false_starts) > 0 else 0
        else:
            longest_false_seq = 0

        # Number of islands of False entries (bounded by True or at boundaries)
        # This means we count all transitions from True to False.
        # Adjusting the island definition to include boundary cases:
        if not arr.all():
            # A False island starts at -1 or at the beginning of the array
            # and ends at 1 or the end of the array
            num_false_islands = len(false_starts)
        else:
            num_false_islands = 0

        return {
            'missing_frag_from_n': consecutive_false_start,
            'missing_frag_from_c': consecutive_false_end,
            'missing_frag_sites': num_false_entries,
            'missing_frag_pct': percentage_false,
            'missing_frag_longest_sequence': longest_false_seq,
            'missing_frag_count': num_false_islands
        }