from typing import List, Optional, Union
import numpy as np
from .psm import PSM

class Spectrum:
    def __init__(self, spectrum_id, **properties):
        self.spectrum_id = spectrum_id
        self.properties = properties
        self.psm_gt: Optional[List[PSM]] = []
        self.psm_candidates: Optional[List[PSM]] = []  # List to hold multiple PSMs associated with this spectrum

    def __repr__(self):
        str_repr = []
        str_repr.append(f'Spectrum ID: {self.spectrum_id}')
        str_repr.append('Ground-truths:')
        str_repr.append('--------------')
        for i, psm in enumerate(self.psm_gt):
            str_repr.append('\t1. {} ({})'.format(
                self.psm.peptide_evidence,
                self.psm.scores
            ))

        str_repr.append('Candidates:')
        str_repr.append('-----------')

        for i, psm_candidate in enumerate(self.psm_candidates):
            str_repr.append("\n\t{}. {} ({})".format(
                i,
                psm_candidate.peptide_evidence,
                psm_candidate.scores
            ))
        return str_repr
    
    def __len__(self):
        return len(self.psm_candidates)

    def add_psm(self, psm: PSM, is_ground_truth=False):
        if is_ground_truth:
            self.psm_gt.append(psm)
        else:
            self.psm_candidates.append(psm)
    
    def rerank(self, score_name, engines, rank_metadata):

        psms = []
        scores = []

        for psm in self.psm_candidates:
            if psm.engine_name in engines:
                scores.append(-psm.scores.get_score(score_name))
                psms.append(psm)
        
        order_idx = np.argsort(scores)

        for i, order_id in enumerate(order_idx):
            psms[order_id].metadata[rank_metadata] = i

    def get_psms_by_engine(self, engine_name):
        return [psm for psm in self.psm_candidates if psm.engine_name == engine_name]
    
    def compare_gt(self, metadata_score, refinements=None, ignore_score=False):

        # TODO: fix
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