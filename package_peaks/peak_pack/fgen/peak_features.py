from ms2rescore.feature_generators.base import FeatureGeneratorBase
from typing import List
from psm_utils import PSM, PSMList

class PeakFeatures(FeatureGeneratorBase):
    """MS2Rescore type feature generator for peak-type features."""

    def __init__(self, *args, **kwargs):
        """Initialize feature generator class."""
        super().__init__(*args, **kwargs)

    @property
    def feature_names(self) -> List[str]:
        """Names of features added to rescoring_features dict."""
        return [
            "peak_tic",
            "peak_count"
        ]
    
    @property
    def name(self) -> str:
        return "general"
    @property
    def input_type(self) -> str:
        return "psm_list"
    
    def add_features(self, psm_list: PSMList, spectra: List) -> None:
        """Compute and add rescoring features to psmlist."""

        for psm, spectrum in zip(psm_list, spectra):
            intensity_arr = spectrum["intensity array"]
            psm.rescoring_features.update({
                "tic": sum(intensity_arr),
                "peak_count": len(intensity_arr),
            })