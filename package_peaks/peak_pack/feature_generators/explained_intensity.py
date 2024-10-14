from ms2rescore.feature_generators.base import FeatureGeneratorBase
from typing import List
from psm_utils import PSM, PSMList

from ..utils import ion_dict_to_matrix, mask_duplicates
from ..annotation.annotation import SpectrumVector

import numpy as np

class ExplainedIntensityFeatures(FeatureGeneratorBase):
    """MS2Rescore type feature generator for peak-type features."""

    def __init__(self, *args, **kwargs):
        """Initialize feature generator class."""
        super().__init__(*args, **kwargs)

    @property
    def feature_names(self) -> List[str]:
        """Names of features added to rescoring_features dict."""
        return [
            "explained_y_pct",
            "explained_b_pct",
            "explained_by_pct",
            "explained_nl_pct",
            "explained_all_pct",
        ]

    @property
    def name(self) -> str:
        return "explained"
    @property
    def input_type(self) -> str:
        return "spectrum_vector"

    def calculate_features(self, sv) -> dict:
        m_by = ion_dict_to_matrix(
            ion_dict=sv.theoretical_mz[""],
            ion_types="y1 b1 y2 b2".split(),
            n=sv.n
        )
        m_all = np.concatenate(
            [
                ion_dict_to_matrix(
                    ion_dict=sv.theoretical_mz[nl],
                    ion_types="y1 b1 y2 b2 a1 a2 c1 c2 x1 x2 z1 z2 p1 p2 p3".split(),
                    n=sv.n
                ) 
                for nl in ["", "-H2O1"]
            ],
            axis=0
        )

        int_by = ion_dict_to_matrix(
            ion_dict=sv.experimental_intensity[""],
            ion_types="y1 b1 y2 b2".split(),
            n=sv.n
        )
        int_all = np.concatenate(
            [
                ion_dict_to_matrix(
                    ion_dict=sv.experimental_intensity[nl],
                    ion_types="y1 b1 y2 b2 a1 a2 c1 c2 x1 x2 z1 z2 p1 p2 p3".split(),
                    n=sv.n
                ) 
                for nl in ["", "-H2O1"]
            ],
            axis=0
        )

        _, mask_by = mask_duplicates(
            m_by,
            preference_list=list(range(4))
        )
        _, mask_all = mask_duplicates(
            m_all,
            preference_list=list(range(4))+list(range(15,19))+list(range(4,15))+list(range(19,30))
        )

        sum_by = np.nansum(np.where(mask_by, int_by, np.nan))
        sum_all = np.nansum(np.where(mask_all, int_all, np.nan))
        sum_nl = sum_all - sum_by
        sum_y, sum_b = int_by.sum(axis=1).reshape(2,2).sum(axis=0)
        arr = np.array(
            [sum_y, sum_b, sum_by, sum_nl, sum_all]
        ) / sv.tic
        features = {k:v for k, v in zip(self.feature_names, arr)}
        return features

    def add_features(self, psm_list: PSMList, spectrum_vector_list: List[SpectrumVector]) -> None:
        """Compute and add rescoring features to psmlist."""
        for psm, spectrum_vector in zip(psm_list, spectrum_vector_list):
            psm.rescoring_features.update(
                self.calculate_features(spectrum_vector)
            )