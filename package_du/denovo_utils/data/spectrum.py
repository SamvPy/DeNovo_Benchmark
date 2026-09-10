from typing import List, Optional, Union
import numpy as np
from .psm import PSM


class Spectrum:
    def __init__(self, spectrum_id, **properties):
        self.spectrum_id = spectrum_id
        self.properties = properties
        self.psm_gold_standard: Optional[List[PSM]] = []
        self.psm_candidates: Optional[List[PSM]] = []  # List to hold multiple PSMs associated with this spectrum

    def __repr__(self):
        str_repr = []
        str_repr.append(f'Spectrum ID: {self.spectrum_id}')
        str_repr.append('Ground-truths:')
        str_repr.append('--------------')
        for i, psm in enumerate(self.psm_gold_standard):
            str_repr.append('\t{}\t{}. {} ({})'.format(
                psm.engine_name,
                i,
                psm.peptide_evidence,
                psm.scores
            ))

        str_repr.append('\nCandidates:')
        str_repr.append('-----------')

        for i, psm_candidate in enumerate(self.psm_candidates):
            str_repr.append("\t{}\t{}. {} ({})".format(
                psm_candidate.engine_name,
                i,
                psm_candidate.peptide_evidence,
                psm_candidate.scores
            ))
        return "\n".join(str_repr)
    
    def __len__(self):
        return len(self.psm_candidates)

    def add_psm(self, psm: PSM, is_gold_standard=False):
        if is_gold_standard:
            self.psm_gold_standard.append(psm)
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
    
    def compare_gold_standard(self, metadata_score, refinements=None, ignore_score=False, ignore_IL=True):

        # TODO: fix
        if self.psm_gold_standard is []:
            return

        for psm in self.psm_candidates:
            for psm_gs in self.psm_gold_standard:
                psm.compare(
                    psm_gt=psm_gs,
                    metadata_score=metadata_score,
                    refinements=refinements,
                    ignore_score=ignore_score,
                    ignore_IL=ignore_IL
                )

    @property
    def engines(self):

        return set([psm.engine_name for psm in self.psm_candidates] + [self.psm_gt.engine_name]) 
