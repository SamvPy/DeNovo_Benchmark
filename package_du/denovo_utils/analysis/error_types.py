from psm_utils import Peptidoform
from peak_pack.annotation.evidence import PeptideEvidence
from nltk import edit_distance
from Levenshtein import editops # To track actions for transpositions
from pyteomics.mass import calculate_mass
import numpy as np
import logging

logger = logging.getLogger(__name__)

MOD_MAPPER = {
    "[+25.980265]-": "b",
    "[Formula:H-2C1O1]-": 'b',
    "[UNIMOD:1]-": "z",
    "[UNIMOD:5]-": "d",
    "[UNIMOD:385]-": "e",
    "Q[UNIMOD:7]": "q",
    "N[UNIMOD:7]": "n",
    "C[UNIMOD:4]": "c",
    "M[UNIMOD:35]": 'v',
    "S[UNIMOD:21]": 's',
    "T[UNIMOD:21]": 't',
    "Y[UNIMOD:21]": 'y',
    "L": "I"
}

MASS_N_TERMINUS = Peptidoform('G').theoretical_mass - (Peptidoform('GGG').theoretical_mass - Peptidoform('GG').theoretical_mass)

def get_mass_tag(tag):
    if tag == '':
        return 0
    elif tag.endswith('-'):
        return Peptidoform(tag+'G/2').theoretical_mass-Peptidoform('G/2').theoretical_mass
    else:
        return Peptidoform(tag+"/2").theoretical_mass - MASS_N_TERMINUS

def mass_difference(s1, s2):
    return get_mass_tag(s1) - get_mass_tag(s2)

def equal_mass(mass_1, mass_2, precision):
    equal = np.round(mass_1 - mass_2, precision) == 0
    return equal

def reverse_opcodes(opcodes):
    """
    Reverse edit opcodes (tag, i, j) that transform s1 → s2 into opcodes for s2 → s1.
    """
    reversed_ops = []

    for tag, i, j in opcodes:
        if tag == 'insert':
            # To undo an insert at s1[i] using s2[j], we delete from s2 at position j
            reversed_ops.append(('delete', j, i))
        elif tag == 'delete':
            # To undo a delete at s1[i], we insert into s2 at position j (which corresponds to s1[i])
            reversed_ops.append(('insert', j, i))
        elif tag == 'replace':
            # To undo a replace, just do it the other way around
            reversed_ops.append(('replace', j, i))
        else:
            raise ValueError(f"Unknown tag: {tag}")

    return reversed_ops

class SequenceComparison:
    def __init__(
            self,
            identifier, 
            peptide_source: Peptidoform,
            peptide_target: Peptidoform,
            mod_mapper=MOD_MAPPER
        ):
        """
        The sequences comparisons are collected in the self.collection object.

        A more pretty format is stored in the self.annotations object

        Explanation of the dictionaries:
        `forward`: when changing the target sequence from the source (source -> target).
        Within the `forward` values, `target` refers to how the target should look like (so what originally was the source).
        `source` refers to how the target used to look like.
        """
        self.identifier = identifier
        self.mod_mapper = mod_mapper
        self.mod_mapper_reverse = {
            v: k for k, v in mod_mapper.items()
        }
        _ = self.mod_mapper_reverse.pop('I')
        self.peptidoform_source = peptide_source 
        self.peptidoform_target = peptide_target
        self.pepstring_source = self._peptidoform_to_lstring(peptide_source, mod_mapper)
        self.pepstring_target = self._peptidoform_to_lstring(peptide_target, mod_mapper)
        self.evidence_added = False
        self.distances = {}
        self.collection = {}
        self.annotations = []
        
    def __repr__(self):
        return f"{self.peptidoform_source.modified_sequence}\n" \
        f"{self.peptidoform_target.modified_sequence}\n" \
        f"{self.annotations}"

    def _add_locations(self, collection, peptide_evidence, opcodes):
        # get ambiguous sites
        ambiguous_indices = []
        ambiguous_tags = peptide_evidence.get_ambiguous_tag_idx()
        for idx_1, idx_2 in ambiguous_tags:
            ambiguous_indices += list(range(idx_1, idx_2+1))
        
        for change in collection:
            ambiguous_site = True # True if all fragment ions are missing
            all_annotated = True # True if all fragment ions are present
            for idx in change['indices']:
                if idx not in ambiguous_indices:
                    ambiguous_site = False
                if idx in ambiguous_indices:
                    all_annotated = False

            # If all fragments are annotated
            if all_annotated:
                change['annotation'] = 'complete'
            elif ambiguous_site:
                change['annotation'] = 'absent'
            else:
                change['annotation'] = 'partly'
            

            change['n_term'] = change['indices'][0]
            change['c_term'] = len(peptide_evidence.peptidoform) - change['indices'][-1] - 1

            if peptide_evidence.peptidoform.properties['n_term'] is not None:
                change['c_term'] += len(peptide_evidence.peptidoform.properties['n_term'])

            action = None
            for opcode in opcodes:
                if opcode[-1] == change['indices'][-1]:
                    action = opcode[0]
            if action == 'delete':
                change['c_term'] += 1
        
        return collection

    def _pretty_parse_tags(self, collections, with_evidence):
        annotations = {'forward': [], 'backward': []}
        for peptide, collection in collections.items():
            for change in collection:
                target = self._peptidoform_to_lstring(change['target'], self.mod_mapper_reverse)
                source = self._peptidoform_to_lstring(change['source'], self.mod_mapper_reverse)
                
                if with_evidence and self.symmetric:
                    annotations[peptide].append(
                        (f'{source} -> {target}', change['ambiguous'])
                    )
                else:
                    annotations[peptide].append(
                        (f'{source} -> {target}', None)
                    )

                md = get_mass_tag(target) - get_mass_tag(source)
                change['mass_difference'] = md
        return annotations

    def _distance_is_symetric(self, c1, c2):
        if len(c1) != len(c2):
            # print(f"WARNING!\n{c1}\n{c2}")
            logger.warning(f'{self.identifier}: Found unequal length tag collection!')
            return False
        
        for (i1, i2) in zip(c1, c2):
            if i1['target'] != i2['source']:
                # print(f"WARNING!\n{c1}\n{c2}")
                logger.warning(f'{self.identifier}: Found asymetric tag collection!')
                return False
        return True

    def _peptidoform_to_lstring(self, peptidoform: Peptidoform, mapper: dict):
        """Parse peptidoforms so modifications are encoded by a single token. (depends on mapper)"""
        if isinstance(peptidoform, Peptidoform):
            new_peptidoform = peptidoform.proforma
        else:
            new_peptidoform = peptidoform
        for p, l in mapper.items():
            new_peptidoform = new_peptidoform.replace(p, l)
        return new_peptidoform.split('/')[0]
    
    def _infer_ambiguity(self, collection_forward, collection_backward):
        """
        Evaluate whether the differing sites are ambiguous between the two sequences.
        
        Can only be run after annotating the tags with `add_locations` with peptide evidence.
        """
        for (tag_1, tag_2) in zip(collection_forward, collection_backward):
            # If both are fully annotated or both have no annotations
            if tag_1['annotation']=='partly' or tag_2['annotation']=='partly':
                tag_1['ambiguous'] = False
                tag_2['ambiguous'] = False
            elif tag_1['annotation']==tag_2['annotation']:
                tag_1['ambiguous'] = True
                tag_2['ambiguous'] = True
            else:
                tag_1['ambiguous'] = False
                tag_2['ambiguous'] = False

    def to_dict(self):
        return {
            'index': self.identifier,
            'levenshtein': self.distances['levenshtein'],
            'damerau-levenshtein': self.distances['damerau-levenshtein'],
            'DL-norm': self.distances['damerau-levenshtein']/max(len(self.pepstring_source), len(self.pepstring_target)),
            'L-norm': self.distances['levenshtein']/max(len(self.pepstring_source), len(self.pepstring_target)),
            'mass_difference': self.distances['mass'],
            'n_sites': self.sites,
            'n_sites_merged': self.sites_merged,
            'sites_pruned': self.sites_pruned,
            'collection': self.collection,
            'annotations': self.annotations,
            'collection_merged': self.collection_merged,
            'annotations_merged': self.annotations_merged,
            'symmetric': self.symmetric
        }

    def add_evidence(self, evidence_source: PeptideEvidence, evidence_target: PeptideEvidence):
        """Add peptide evidence."""
        self.pe_source = evidence_source
        self.pe_target = evidence_target

    def calculate_distance(self):
        """Calculate the (Damerau-)Levenshtein distance. Required before annotating differences."""
        self.distances['damerau-levenshtein'] = edit_distance(
            s1=self.pepstring_source,
            s2=self.pepstring_target,
            transpositions=True
        )
        self.opcodes_forward = editops(
            self.pepstring_source,
            self.pepstring_target
        )
        self.opcodes_backward = reverse_opcodes(
            self.opcodes_forward
        )
        self.distances['levenshtein'] = len(self.opcodes_backward)
        self.distances['mass'] = self.peptidoform_target.theoretical_mass - self.peptidoform_source.theoretical_mass

    def annotate_tags(self, with_evidence=False):
        """
        Parse the differences between the two sites in an easily readible fashion.
        
        Also evaluates whether the differences are ambiguous if with_evidence is True.
        """
        assert len(self.distances) == 3

        # 1. Parse the opcodes towards tags
        tracker_1 = ChangeTracker()
        tracker_2 = ChangeTracker()
        collection_forward = tracker_1.track_changes(
            s1=self.pepstring_source,
            s2=self.pepstring_target,
            opcodes=self.opcodes_forward
        )
        collection_backward = tracker_2.track_changes(
            s1=self.pepstring_target,
            s2=self.pepstring_source,
            opcodes=self.opcodes_backward
        )

        # 2. If there is spectral evidence, annotate the following
        # - mass differences
        # - Distance from terminus
        # - Spectral evidence in both tags present or absent (ambiguity)
        if with_evidence:
            collection_forward = self._add_locations(
                collection=collection_forward,
                peptide_evidence=self.pe_target,
                opcodes=self.opcodes_forward
            )
            collection_backward = self._add_locations(
                collection=collection_backward,
                peptide_evidence=self.pe_source,
                opcodes=self.opcodes_backward
            )
        # 3. Collect the collections
        self.collection = {'forward': collection_forward, 'backward': collection_backward}

        # 4. Check if everything went as expected
        symmetric = self._distance_is_symetric(c1=collection_forward, c2=collection_backward)
        self.symmetric = symmetric
        
        # Sometimes, this is not symmetric...
        # If not symmetric, ambiguity cannot be defined. Should never occur and will produce a warning.
        # 5. Infer ambiguity based on peak annotations.
        if symmetric:
            if with_evidence:
                self._infer_ambiguity(collection_forward=collection_forward, collection_backward=collection_backward)

        # 5. Make the tags and differences print-ready
        # This _pretty_parse_tags also adds mass difference annotations to the tags
        # TODO: Separate this functionality
        self.annotations = self._pretty_parse_tags(collections=self.collection, with_evidence=with_evidence)

        # 6. Try to combine tags within the collections 
        # Requires the self.collection attribute!
        tm = TagMerger()
        collection_merged = tm.merge_tag_workflow(self)
        self.collection_merged = {}
        if with_evidence:
            # Add type of peak annotations, locations from N-C terminus etc.
            self.collection_merged['forward'] = self._add_locations(
                collection=collection_merged['forward'],
                peptide_evidence=self.pe_target,
                opcodes=self.opcodes_forward
            )
            self.collection_merged['backward'] = self._add_locations(
                collection=collection_merged['backward'],
                peptide_evidence=self.pe_source,
                opcodes=self.opcodes_backward
            )
            # Based on peak annotation, are the tags ambiguous?
            self._infer_ambiguity(
                collection_forward=self.collection_merged['forward'],
                collection_backward=self.collection_merged['backward']
            )

        # 7. Parse the collected tags into a nicely printable format
        # This also adds mass differences for the tags
        self.annotations_merged = self._pretty_parse_tags(collections=self.collection_merged, with_evidence=with_evidence)

        # 8. Count number of tags
        self.sites = len(self.annotations['forward'])
        self.sites_merged = len(self.annotations_merged['forward'])
        self.sites_pruned = self.sites - self.sites_merged

class ChangeTracker:
    def __init__(self):
        self.change_s1 = ''
        self.change_s2 = ''
        self.collection = []
        self.indices = []
        self.prev_action = None
    

    def _collect_change(self):
        '''Collect an entry.'''
        self.collection.append(
            {
                'source': self.change_s1,
                'target': self.change_s2,
                'indices': self.indices
            }
        )
        self._flush()
    
    def _flush(self):
        '''Clear the tracker.'''
        self.change_s1 = ''
        self.change_s2 = ''
        self.indices = []
        self.prev_action = None

    def track_changes(self, s1, s2, opcodes, **kwargs):
        self.collection = []

        deleted = 0
        prev_index = 0
        flag_start = True

        for (action, idx_s1, idx_s2) in opcodes:

            # # New site of change
            # if action == 'delete':
            #     deleted += 1
            # idx_s2_deletion = idx_s2 - deleted

            # if idx_s2_deletion - prev_index > 1 and not flag_start:
            #     # If multiple deletes after each other
            #     if len(self.indices)>0 and (self.indices[-1]==idx_s2):
            #         pass
            #     else:
            #         print('passed here')
            #         self._collect_change()

            # Detect whether a new site starts
            if flag_start or len(self.indices) == 0:
                pass
            elif action == 'delete':
                # Multiple deletes after each other
                if (idx_s2 == prev_index):
                    pass
                # Transposition
                elif (self.prev_action == 'insert') and (idx_s2 == prev_index+2):
                    pass
                # Transposition case 2
                elif (self.prev_action == 'replace') and (idx_s2 == prev_index+1):
                    pass
                # Otherwise, this is a new tag
                else:
                    self._collect_change()
            else:
                if idx_s2 - prev_index > 1:
                    self._collect_change()


            if action == 'insert':
                # change_s1 += s1[idx_s1]
                if self.prev_action == 'delete':
                    self.change_s1 += s1[idx_s1-1]
                    self.change_s2 += s2[idx_s2-1]
                    self.change_s2 += s2[idx_s2]

                else:
                    self.change_s2 += s2[idx_s2]

            elif action == 'replace':
                # # Fixes weird bug where it does not detect new site
                if self.prev_action == 'delete':
                    if self.indices[-1] != idx_s2:
                        self._collect_change()
                self.change_s1 += s1[idx_s1]
                self.change_s2 += s2[idx_s2]

            elif action == 'delete':
                # Transposition
                if self.prev_action == 'insert':
                    # print(idx_s1, idx_s2, idx_s2_deletion, deleted)
                    self.change_s1 += s1[idx_s1-1]
                    self.change_s2 += s2[idx_s2-1]
                    self.change_s1 += s1[idx_s1]
                
                elif self.prev_action == 'replace':
                    # print(idx_s1, idx_s2, idx_s2_deletion, prev_index, self.indices)
                    if self.indices[-1]+1 != idx_s2:
                        self._collect_change()
                    self.change_s1 += s1[idx_s1]

                elif self.prev_action == 'delete' and not (self.indices[-1]==idx_s2):   
                    self.change_s1 += s1[idx_s1]

                else:
                    self.change_s1 += s1[idx_s1]
                # change_s2 += s2[idx_s2-deleted]
            

            if action == 'delete' and self.prev_action == 'insert':  
                self.indices.append(idx_s2-1)
            else:
                self.indices.append(idx_s2)

            # prev_index = idx_s2_deletion
            prev_index = idx_s2
            self.prev_action = action
            flag_start = False

        self._collect_change()
        return self.collection


class TagMerger:
    def __init__(self):
        """
        Combines tags with mass differences.

        Some tags inferred by the ChangeTracker are not a single tag.
        They are separated by a correct subsequence which happened to be correct.
        Nonetheless, the mass shift introduced by the mismatch (tag), shifts the location
        of the peaks in the mass spectrum. Therefore, this function merges tags separated
        by a correct subsequence, which together, forms a mass shift of 0 (isobaric) or
        equal to the precursor mass difference if there is one.

        Methods
        -------
        merge_tags_workflow:
            Main method to run the tag combiner. Requires a SequenceComparison object.
        """
        pass

    def _final_index_delete(self, opcodes, final_idx):
        final_action = None
        for opcode in opcodes:
            if opcode[2] == final_idx:
                final_action = opcode[0]
        return final_action

    def _create_new_collection(self, old_tag_ids, collection, composite_tags):
        '''From the merged tags, create a new collection for the peptides. Also, annotate it again.'''
        
        new_tags = []

        if len(old_tag_ids) == 0:
            return collection

        keep_tag = True
        new_tags_n = len(old_tag_ids)-1
        current_i = 0
        no_more_composites = False

        for tag_i, tag in enumerate(collection):

            if no_more_composites:
                pass
            elif tag_i == old_tag_ids[current_i][0]:
                keep_tag = False

            if keep_tag:
                new_tags.append(tag)
            else:
                if tag_i == old_tag_ids[current_i][1]:
                    new_tags.append(composite_tags[current_i])
                    keep_tag = True
                    current_i += 1
                    if current_i > new_tags_n:
                        no_more_composites = True

        return new_tags
        
    def _construct_tag_list(self, result):
        new_composite_tags_forward = []
        new_composite_tags_backward = []
        for tag_forward, tag_backward, index_forward, index_backward in zip(
            result['forward']['tags'],
            result['backward']['tags'],
            result['forward']['indices'],
            result['backward']['indices']
        ):
            new_composite_tags_forward.append(
                {
                    'target': tag_backward,
                    'source': tag_forward,
                    'indices': list(range(index_forward[0], index_forward[1]+1))
                }
            )
            new_composite_tags_backward.append(
                {
                    'target': tag_forward,
                    'source': tag_backward,
                    'indices': list(range(index_backward[0], index_backward[1]+1))
                }
            )
        return new_composite_tags_forward, new_composite_tags_backward

    def _prune_amino_acids(self, ids_forward, tags_forward, ids_backward, tags_backward):
        assert len(ids_forward) == len(ids_backward)
        
        ids_forward_new = []
        ids_backward_new = []
        tags_forward_new = []
        tags_backward_new = []

        for id_forward, tag_forward, id_backward, tag_backward in zip(
            ids_forward, 
            tags_forward,
            ids_backward,
            tags_backward
        ):
            
            # Evaluate equal beginning amino acids and prune
            equal_start_aas = 0
            for i in range(min(len(tag_forward), len(tag_backward))):
                if tag_forward[i] == tag_backward[i]:
                    equal_start_aas += 1
                else:
                    break

            # Shift the start of the tag
            id_forward[0] += equal_start_aas
            id_backward[0] += equal_start_aas
            tag_forward = tag_forward[equal_start_aas:]
            tag_backward = tag_backward[equal_start_aas:]

            # Evaluate equal ending amino acids and prune
            equal_end_aas = 0
            for i in range(1, min(len(tag_forward), len(tag_backward))+1):
                if tag_forward[-i] == tag_backward[-i]:
                    equal_end_aas += 1
                else:
                    break
            
            # Shift the end of the tag
            id_forward[1] -= equal_end_aas
            id_backward[1] -= equal_end_aas
            if equal_end_aas > 0:
                tag_forward = tag_forward[:-equal_end_aas]
                tag_backward = tag_backward[:-equal_end_aas]

            ids_forward_new.append(id_forward)
            ids_backward_new.append(id_backward)
            tags_forward_new.append(tag_forward)
            tags_backward_new.append(tag_backward)

        return [(ids_forward_new, tags_forward_new), (ids_backward_new, tags_backward_new)]

    def _merge_tags(self, comparison: SequenceComparison, peptide_id = 'forward'):
        """
        Infer whether two tags are actually a single tag.

        Sometimes tags can have a mass difference, which becomes 0 or equal to the peptide mass difference when combined.
        Essentially, this is a long tag, which was seperated by a middle part, which happens to be equal between the sequences.

        When this is not done, some cases will show up as being non-symmetric, eventhough the symmetric nature of opcodes is enforced.

        This post-processing step will only be applied to tags which have an actual mass difference, thus who are not isobaric.
        Additionally, these collapsed tags will gain the flag 'Collapsed', in order to prevent extremely long tags.
        """
        searching_for_merge = False
        current_difference = 0
        ignore_current_mass_difference = False

        new_tag_idxs_list = []
        old_tag_idxs_list = []

        if peptide_id == 'backward':
            precursor_md = -np.round(comparison.distances['mass'], 2)
            pepstring = comparison.pepstring_source
            opcodes = comparison.opcodes_backward
        else:
            precursor_md = np.round(comparison.distances['mass'], 2)
            pepstring = comparison.pepstring_target
            opcodes = comparison.opcodes_forward
        
        precursor_md_found = equal_mass(precursor_md, 0, 2)

        for tag_i, tag in enumerate(comparison.collection[peptide_id]):
            # First, evaluate whether this tag has the same mass as the precursor difference
            # If so, we will see this as a tag which is isobaric
            # In reality, this will make no sense, to do it this way, but I want to extract these errors specifically.
            if not precursor_md_found:

                # If we find a tag with the same mass difference as the precursor mass difference
                if equal_mass(tag['mass_difference'], precursor_md, 2):

                    # No longer search for this kind of mass difference
                    precursor_md_found = True
                    # And allow the current cycle to ignore this mass difference
                    # while not breaking a potential larger tag chain
                    ignore_current_mass_difference = True

                # Else we could identify a combined mass of several tags as equal the precursor mass difference.
                # This we will see as finding a new tag!
                elif equal_mass(tag['mass_difference']+current_difference, precursor_md, 2):
                    # We will only arrive here when we are searching for a merge (as current difference will be !=0)
                    end_tag = tag['indices'][-1]
                    
                    # if self._final_index_delete(opcodes=opcodes, final_idx=end_tag) != 'delete':
                    #     end_tag += 1

                    new_tag_idxs_list.append([start_new_tag, end_tag])
                    old_tag_idxs_list[-1].append(tag_i)

                    # and of course, stop searching
                    searching_for_merge = False
                    current_difference = 0
                    ignore_current_mass_difference = True

            # Isobaric tag: We can ignore this.
            if (np.round(tag['mass_difference'], 2) == 0) or (ignore_current_mass_difference):
                ignore_current_mass_difference = False
                continue

            # A tag with a mass difference! We will try to merge this one or start a new potential tag
            else:
                # If we were searching already, we will check if the mass difference becomes 0 for
                # this larger tag
                if searching_for_merge:
                    current_difference += tag['mass_difference']
                    # Check if we found the end
                    
                    # 1. We found a larger tag!
                    if equal_mass(current_difference, 0, 2):
                        end_tag = tag['indices'][-1]
                        
                        # if self._final_index_delete(opcodes=opcodes, final_idx=end_tag) != 'delete':
                        #     end_tag += 1

                        new_tag_idxs_list.append([start_new_tag, end_tag])
                        old_tag_idxs_list[-1].append(tag_i)

                        # and of course, stop searching
                        searching_for_merge = False
                        current_difference = 0

                    # 2. We did not find the end, so we need to continue searching
                    else:
                        continue
                
                else:
                    # We were not searching already. So this could be the start of a new, larger tag.
                    current_difference += tag['mass_difference']
                    searching_for_merge = True
                    start_new_tag = tag['indices'][0]
                    old_tag_idxs_list.append([tag_i])

        new_tags = []
        for tag_0, tag_1 in new_tag_idxs_list:
            if self._final_index_delete(opcodes=opcodes, final_idx=tag_1) != 'delete':
                new_tags.append(pepstring[tag_0: tag_1+1])
            else:
                if tag_0 == tag_1:
                    new_tags.append(pepstring[tag_0])
                else:
                    new_tags.append(pepstring[tag_0: tag_1])

        # new_tags = [pepstring[tag[0]: tag[1]] for tag in new_tag_idxs_list]

        # If we were still searching for something to merge, don't merge them
        # Don't know why this would occur, but apparently happens...
        if len(old_tag_idxs_list)>0:
            if len(old_tag_idxs_list[-1])==1:
                if len(old_tag_idxs_list)==1:
                    old_tag_idxs_list = []
                else:
                    old_tag_idxs_list = old_tag_idxs_list[:-1]


        return new_tag_idxs_list, new_tags, old_tag_idxs_list
    
    def merge_tag_workflow(self, comparison):

        # Merge tags for forward
        ids_forward, tags_forward, old_tag_ids_forward = self._merge_tags(comparison, peptide_id='forward')

        # Merge tags for backward
        ids_backward, tags_backward, old_tag_ids_backward = self._merge_tags(comparison, peptide_id='backward')

        # Prune equal amino acids
        (ids_forward, tags_forward), (ids_backward, tags_backward) = self._prune_amino_acids(
            ids_forward=ids_forward,
            ids_backward=ids_backward,
            tags_forward=tags_forward,
            tags_backward=tags_backward
        )

        # Construct new tags
        result = {
            'forward': {'indices': ids_forward, 'tags': tags_forward},
            'backward': {'indices': ids_backward, 'tags': tags_backward}
        }
        new_composite_tags_forward , new_composite_tags_backward = self._construct_tag_list(result)

        # Create a new collection with the new tags
        new_tags_forward = self._create_new_collection(
            old_tag_ids_forward,
            comparison.collection['forward'],
            new_composite_tags_forward
        )
        new_tags_backward = self._create_new_collection(
            old_tag_ids_backward,
            comparison.collection['backward'],
            new_composite_tags_backward
        )

        # I need to annotate the sequence comparison here !
        return {'forward': new_tags_forward, 'backward': new_tags_backward}