from typing import List, Optional, Union
import numpy as np
from .psm import PSM

class Spectrum:
    def __init__(self, spectrum_id, **properties):
        self.spectrum_id = spectrum_id
        self.properties = properties
        self.psm_gt: Optional[PSM] = None
        self.psm_candidates: Optional[List[PSM]] = []  # List to hold multiple PSMs associated with this spectrum

    def __repr__(self):
        str_repr = "Spectrum ID: {}\nGround-truth: {} ({})".format(
            self.spectrum_id,
            self.psm_gt.peptide_evidence,
            self.psm_gt.scores
        )
        str_repr += "\nCandidates:"
        for psm_candidate in self.psm_candidates:
            str_repr += "\n\t{} ({})".format(
                psm_candidate.peptide_evidence,
                psm_candidate.scores
            )
        return str_repr
    
    def __len__(self):
        return len(self.psm_candidates)

    def add_psm(self, psm: PSM, is_ground_truth=False):
        if is_ground_truth:
            self.psm_gt = psm
        else:
            self.psm_candidates.append(psm)
    
    def rerank(self, score_name, engines):

        psms = []
        scores = []

        for psm in self.psm_candidates:
            if psm.engine_name in engines:
                scores.append(-psm.scores.get_score(score_name))
                psms.append(psm)
        
        order_idx = np.argsort(scores)

        for i, order_id in enumerate(order_idx):
            psms[order_id].rank = i
            psms[order_id].metadata['previous_rank'] = psm.rank

    def get_psms_by_engine(self, engine_name):
        return [psm for psm in self.psm_candidates if psm.engine_name == engine_name]
    
    def compare_gt(self, metadata_score, refinements=None, ignore_score=False):

        if self.psm_gt is None:
            return

        for psm in self.psm_candidates:
            psm.compare(
                psm_gt=self.psm_gt,
                metadata_score=metadata_score,
                refinements=refinements,
                ignore_score=ignore_score
            )

    @property
    def engines(self):

        return set([psm.engine_name for psm in self.psm_candidates] + [self.psm_gt.engine_name]) 