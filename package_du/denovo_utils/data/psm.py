from ..analysis.metrics import aa_match, convert_peptidoform, mass_diff
import numpy as np

class PSM:
    def __init__(
            self,
            peptidoform,
            score,
            engine_name,
            rank,
            aa_score = None,
            is_ground_truth=False,
            peptide_evidence=None
        ):
        self.peptidoform = peptidoform  # The identified sequence
        self.scores = Score(score, engine_name)  # Confidence score or probability
        self.scores.add_score(aa_score, engine_name, score_type='aa')
        self.engine_name = engine_name  # The tool that generated the PSM
        self.rank = rank
        self.is_ground_truth = is_ground_truth  # Whether this is the ground truth PSM
        self.peptide_evidence = peptide_evidence
        self.refinement = {}
        self.evaluation = {}
        self.metadata = {}

    def __eq__(self, psm: 'PSM'):
        """Exact string matching of peptidoform."""
        return self.peptidoform.proforma == psm.peptidoform.proforma
    
    def __repr__(self):
        return str(self.to_dict())
    
    @property
    def modifications(self):
        mod_list = []
        if not self.peptidoform.is_modified:
            return mod_list
        
        n_mod = self.peptidoform.properties['n_term']
        if isinstance(n_mod, list):
            mod_list += n_mod

        for aa, mods in self.peptidoform.parsed_sequence:
            if mods:
                mod_list += mods
        return list(set(mod_list))

    def has_modification(self, modification):
        return modification in self.modifications

    def to_dict(self):
        return {
            'peptidoform': self.peptidoform,
            'scores': self.scores,
            'engine_name': self.engine_name,
            'peptide_evidence': self.peptide_evidence,
            'refinement': self.refinement,
            'evaluation': self.evaluation
        }

    def add_peptide_evidence(self, peptide_evidence):
        self.peptide_evidence = peptide_evidence

    def compare(self, psm_gt: 'PSM', metadata_score, refinements=None, tolerance=.02, ignore_score=False):
        # Recursively compare the refinements of the base PSM with ground truth
        if refinements is not None:
            for refinement in refinements:
                if refinement not in self.refinement.keys():
                    continue
                psm = self.refinement[refinement][0]
                if psm is not None:
                    psm.compare(psm_gt, metadata_score, ignore_score=ignore_score)

        # Extract the scores associated with the ground-truth and current PSM
        score = self.scores.get_score(
            metadata=metadata_score,
            ignore_missing_score=ignore_score
        )
        score_gt = psm_gt.scores.get_score(
            metadata=metadata_score,
            ignore_missing_score=ignore_score
        )
        # Perform the evaluation against a ground-truth
        evaluation = Evaluation(gt_name=psm_gt.engine_name, other_name=self.engine_name)
        evaluation.evaluate(psm_gt, self, score_gt, score)

        self.evaluation[metadata_score] = evaluation

    def add_refinement(self, psm: 'PSM', overwrite=False):
        equal_sequence = psm == self
        metadata = psm.engine_name

        if equal_sequence:
            self.scores.add_score(
                score=psm.scores.get_score(psm.engine_name),
                metadata=metadata,
                score_type='peptide',
                overwrite=overwrite
            )
            # Prevent saving duplicate psms if postprocessor didnt change the PSM
            psm=None

        self.refinement[metadata] = (psm, equal_sequence)


class Score:
    def __init__(self, score=None, metadata=None, score_type='peptide'):
        self.score_dict = {'peptide': {}, 'aa': {}}
        if score is not None and metadata is not None:
            self.add_score(score, metadata, score_type)
        
    def __repr__(self):
        return str(self.score_dict)
    
    def add_score(self, score, metadata, score_type, overwrite=False):
        if score is None:
            return

        if (not overwrite) and \
            (metadata in self.score_dict[score_type].keys()) and \
            (self.get_score(metadata, score_type)!=score):
            raise Exception(f'{metadata} score is already registered as {self.score_dict[score_type][metadata]}.\
                            Set overwrite to True to overwrite.')
        
        if isinstance(score, str):
            score = eval(score)
        self.score_dict[score_type][metadata] = score

    def get_score(self, metadata, score_type='peptide', ignore_missing_score=False):
        
        if metadata not in self.score_dict[score_type].keys():
            if ignore_missing_score:
                return None
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
    def __init__(self, gt_name, other_name):
        self.engines = {"gt": gt_name, "other": other_name}
        self.error_type = None
        self.isobaric_aa_errors = []
        self.isobaric_tag_errors = []
        self.aa_match = None
        self.score_diff = None

    def __repr__(self):
        return str((self.error_type, self.score_diff))

    def to_dict(self):
        return {
            "error_type": self.error_type,
            "isobaric_aa_errors": self.isobaric_aa_errors,
            "isobaric_tag_errors": self.isobaric_tag_errors,
            "aa_match": self.aa_match,
            "score_diff": self.score_diff
        }

    def evaluate(self, psm_1, psm_2, score_1, score_2):
        evaluate_isobars = True

        if isinstance(score_1, float) and isinstance(score_2, float):
            self.score_diff = round(score_2 - score_1, 3)

        if psm_1.peptidoform == psm_2.peptidoform:
            self.error_type = 'match'
            return

        _, pep_match, (match_1, match_2), iso_errs, tols = aa_match(
            convert_peptidoform(psm_1.peptidoform),
            convert_peptidoform(psm_2.peptidoform),
            aa_dict={}
        )
        self.isobaric_aa_errors = iso_errs

        if pep_match:
            self.error_type = 'isobaric_aa'
            return

        pe_1 = psm_1.peptide_evidence
        pe_2 = psm_2.peptide_evidence
        if pe_1 is None:
            evaluate_isobars = False

        if pe_2 is None:
            evaluate_isobars = False

        precursor_mass_diff = abs(mass_diff(psm_1.peptidoform.theoretical_mass, psm_2.peptidoform.theoretical_mass, True)) >= 0.5

        if evaluate_isobars and not precursor_mass_diff:
            match_evidence_1 = evaluate_isobaricity(pe_1, match_1)
            match_evidence_2 = evaluate_isobaricity(pe_2, match_2)

            if match_evidence_1 and match_evidence_2:
                tag_list_1 = pe_1.ambiguous_tags
                tag_list_2 = pe_2.ambiguous_tags

                try:
                    assert len(tag_list_1) == len(tag_list_2)
                except AssertionError:
                # This only happens when the tolerance is a little bit of
                # Due to the difference in tolerance calculation in aa_match and rustyms

                # In aa_match: Tolerance is calculated between the two hits
                # In rustyms/peptide-evidence: Tolerance calculated between observed and a given hit
                # Thus, if peak of hit 1 has tolerance with observed of 19 and hit 2 21, these two will
                # be considered equal by aa_match, yet have different evidence according to peptide-evidence.
                # Meaning these seemingly equal sequences will be differently scored as well as the PSM-features
                # will be considerably different.

                # In this special case, which seems to happen 1 in 5.000, consider it a non-isobaric hit
                    # Ugly duplicate code
                    if self.score_diff is None:
                        self.error_type = "Different"
                    elif self.score_diff < 0:
                        self.error_type = 'Lower score'
                    else:
                        self.error_type = 'Higher score'
                    return

                # issue with the tag creator ?
                for tag_1, tag_2 in zip(tag_list_1, tag_list_2):
                    self.isobaric_tag_errors.append(f'{tag_1}/{tag_2}')
                    self.error_type = 'isobaric_peak'
                
                self.error_type = "isobaric_peak"
                return

        if self.score_diff is None:
            self.error_type = "Different"
        elif self.score_diff < 0:
            self.error_type = 'Lower score'
        else:
            self.error_type = 'Higher score'
        
        if self.error_type is None:
            raise Exception(f"Error type is not properly set")

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