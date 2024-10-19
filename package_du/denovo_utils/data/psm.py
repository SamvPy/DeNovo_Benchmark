from ..analysis.metrics import aa_match, convert_peptidoform
import numpy as np

class PSM:
    def __init__(
            self,
            peptidoform,
            score,
            engine_name,
            aa_score = None,
            is_ground_truth=False,
            peptide_evidence=None
        ):
        self.peptidoform = peptidoform  # The identified sequence
        self.scores = Score(score, engine_name)  # Confidence score or probability
        self.scores.add_score(aa_score, engine_name, score_type='aa')
        self.engine_name = engine_name  # The tool that generated the PSM
        self.is_ground_truth = is_ground_truth  # Whether this is the ground truth PSM
        self.peptide_evidence = peptide_evidence
        self.refinement = {}
        self.evaluation = {}

    def __eq__(self, psm: 'PSM'):
        """Exact string matching of peptidoform."""
        return self.peptidoform.proforma == self.psm.proforma

    def add_peptide_evidence(self, peptide_evidence):
        self.peptide_evidence = peptide_evidence

    def compare(self, psm_gt: 'PSM', metadata_score, refinement=None, tolerance=.02):
        if refinement is not None:
            psm = self.refinement[refinement][0]
            psm.compare(psm_gt, metadata_score)

        score = self.scores.get_score(metadata=metadata_score)
        score_gt = psm_gt.scores.get_score(metadata=metadata_score)
        evaluation = Evaluation()
        evaluation.evaluate(psm_gt, self, score_gt, score)
        self.evaluation[metadata_score] = evaluation

    def add_refinement(self, psm: 'PSM', metadata):
        equal_sequence = psm == self
        if equal_sequence:
            self.scores.add_score(
                score=psm.scores.get_score(psm.engine_name),
                metadata=psm.engine_name,
                score_type='peptide'
            )

        self.refinement[metadata] = (psm, equal_sequence)


class Score:
    def __init__(self, score=None, metadata=None, score_type='peptide'):
        self.score_dict = {}
        if score is not None and metadata is not None:
            self.add_score(score, metadata, score_type)
    
    def add_score(self, score, metadata, score_type, overwrite=False):
        if score is None:
            return

        if (not overwrite) and \
            (metadata in self.score_dict[score_type].keys()) and \
            (self.get_score(metadata, score_type)!=score):
            raise Exception(f'{metadata} score is already registered as {self.score_dict[score_type][metadata]}.\
                            Set overwrite to True to overwrite.')
        
        self.score_dict[score_type][metadata] = score

    def get_score(self, metadata, score_type='peptide'):
        if metadata not in self.score_dict[score_type].keys():
            raise Exception(f'Only scores from {list(self.score_dict[score_type].keys())} are registered.')
        return self.score_dict[score_type].get(metadata)




class Evaluation:
    """
    Evaluator object comparing two peptidoforms.

    The match type can be the following four values:
    - match: Exact sequence match
    - isobaric_aa: Exact sequence match on mass level (ignoring aa with unresolvable mass)
    - isobaric_peak: Exact sequence match while allowing errors due to incomplete ion coverage
    - mismatch: Sequence prediction does not match
    """
    def __init__(self):
        self.error_type = None
        self.isobaric_aa_errors = []
        self.isobaric_tag_errors = []
        self.aa_match = None
        self.score_diff = None

    def evaluate(self, psm_1, psm_2, score_1, score_2):
        skip_isobars = True
        self.score_diff = score_2 - score_1

        if psm_1.peptidoform == psm_2.peptidoform:
            self.error_type = 'match'

        _, pep_match, (match_1, match_2), iso_errs = aa_match(
            convert_peptidoform(psm_1.peptidoform),
            convert_peptidoform(psm_2.peptidoform)
        )
        self.isobaric_aa_errors = iso_errs

        if pep_match:
            self.error_type = 'isobaric_aa'
            return

        if pe_1 is None:
            pe_1 = psm_1['metadata']['peptide_evidence']
            skip_isobars = False

        if pe_2 is None:
            pe_2 = psm_2['metadata']['peptide_evidence']
            skip_isobars = False

        if not skip_isobars:
            match_evidence_1 = evaluate_isobaricity(pe_1, match_1)
            match_evidence_2 = evaluate_isobaricity(pe_2, match_2)

            if match_evidence_1 and match_evidence_2:
                tag_list_1 = pe_1.ambiguous_tags
                tag_list_2 = pe_2.ambiguous_tags
                assert len(tag_list_1) == len(tag_list_2)
                for tag_1, tag_2 in zip(tag_list_1, tag_list_2):
                    self.isobaric_tag_errors.append(f'{tag_1}/{tag_2}')
                    self.error_type = 'isobaric_peak'

        if self.score_diff < 0:
            self.error_type = 'Lower score'
        else:
            self.error_type = 'Higher score'
        

def evaluate_isobaricity(pe, aa_match_array):

    # Last element in match array is precursor mass matching
    if pe.peptidoform.properties['n_term'] is not None:
        # N-term mods are seperate match index.
        aa_match_array = aa_match_array[1:-1] 
    else:
        aa_match_array = aa_match_array[:-1] 
    
    return aa_match_array[pe.evidence].all()

def filter_out_indices(sequence, index_pairs):
    # Sort the index pairs by start index for efficiency
    index_pairs = sorted(index_pairs)
    
    filtered_sequence = []
    current_position = 0
    
    for start, end in index_pairs:
        # Append the sequence from the current position to the start of the range to be removed
        filtered_sequence.extend(sequence[current_position:start])
        # Move the current position to the end of the range to be removed
        current_position = end + 1  # Assuming end is inclusive
    
    # Append any remaining elements after the last index pair
    filtered_sequence.extend(sequence[current_position:])
    
    return np.array(filtered_sequence)