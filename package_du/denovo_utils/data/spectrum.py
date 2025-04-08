from typing import List, Optional, Union
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

        rank_dict = {}
        for psm in self.psm_candidates:
            if psm.engine_name in engines:
                rank_dict[psm] = psm.scores.get_score(score_name)
        
        rank_dict = {
            psm: score for psm, score in sorted(rank_dict.items(), key=lambda item: item[1])
        }

        for i, psm in enumerate(rank_dict.keys()):
            psm['metadata']['previous_rank'] = psm['rank']
            psm['rank'] = 1

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