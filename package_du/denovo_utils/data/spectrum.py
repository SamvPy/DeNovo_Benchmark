from typing import List, Optional, Union
from .psm import PSM

class Spectrum:
    def __init__(self, spectrum_id):
        self.spectrum_id = spectrum_id
        self.psm_gt: Optional[PSM] = None
        self.psm_candidates: Optional[List[PSM]] = []  # List to hold multiple PSMs associated with this spectrum

    def add_psm(self, psm: PSM, is_ground_truth=False):
        if is_ground_truth:
            self.psm_gt = psm
        else:
            self.psm_candidates.append(psm)
    
    def get_psms_by_engine(self, engine_name):
        return [psm for psm in self.psm_candidates if psm.engine_name == engine_name]
    
    def compare_gt(self, metadata_score, refinement=None):

        if self.psm_gt is None:
            return

        for psm in self.psm_candidates:
            psm.compare(
                psm_gt=self.psm_gt,
                metadata_score=metadata_score,
                refinement=refinement,
            )