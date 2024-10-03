from ms2rescore.feature_generators.base import FeatureGeneratorBase
from typing import List
from psm_utils import PSM, PSMList
from ..utils import ion_dict_to_matrix
from ..annotation.annotation import SpectrumVector
import numpy as np
from tqdm import tqdm


class PPMFeatures(FeatureGeneratorBase):
    """MS2Rescore type feature generator for ppm-errors"""
    def __init__(self, *args, **kwargs):
        """Initialize feature generator class."""
        super().__init__(*args, **kwargs)

    @property
    def feature_names(self) -> List[str]:
        """Names of features added to rescoring_features dict."""
        return [
            "ppm_mean_y",
            "ppm_mean_b",
            "ppm_mean_by"
        ]

    @property
    def name(self) -> str:
        return "explained"
    @property
    def input_type(self) -> str:
        return "spectrum_vector"

    def add_features(self, psm_list: PSMList, spectrum_vector_list: List[SpectrumVector]) -> None:
        """Compute and add rescoring features to psmlist."""
        for psm, spectrum_vector in tqdm(zip(psm_list, spectrum_vector_list)):
            ppm_diff = spectrum_vector.ppm_diff

            ppm_y = ion_dict_to_matrix(
                ion_dict=ppm_diff,
                ion_types=["y1", "y2"],
                n=spectrum_vector.n
            ).flatten()

            ppm_b = ion_dict_to_matrix(
                ion_dict=ppm_diff,
                ion_types=["b1", "b2"],
                n=spectrum_vector.n
            ).flatten()

            ppm_all = ion_dict_to_matrix(
                ion_dict=ppm_diff,
                ion_types=["y1", "y2", "b1", "b2"],
                n=spectrum_vector.n
            ).flatten()
            
            psm.rescoring_features.update({
                "ppm_mean_y": -1 if np.isnan(ppm_y).all() else np.abs(np.nanmean(ppm_y)),
                "ppm_mean_b": -1 if np.isnan(ppm_b).all() else np.abs(np.nanmean(ppm_b)),
                "ppm_mean_all": -1 if np.isnan(ppm_all).all() else np.abs(np.nanmean(ppm_all)),
                "ppm_precursor": np.abs(spectrum_vector.precursor_ppm)
            })