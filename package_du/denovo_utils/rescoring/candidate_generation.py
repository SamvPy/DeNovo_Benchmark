import random
import numpy as np
import pandas as pd

from itertools import chain, product
from peak_pack.annotation import PeptideEvidence

from pyteomics import mass as ptmass
from pyteomics import proforma

from psm_utils import Peptidoform, PSMList, PSM

class CandidateGenerator:
    def __init__(self, aa_comp, unimod2pyteomics, pyteomics2unimod, tolerance=.02):
        self.tolerance = tolerance
        self.aa_comp = aa_comp
        self.unimod2pyteomics = unimod2pyteomics
        self.pyteomics2unimod = pyteomics2unimod

    def load_mass_table(self, mass_table):
        self.mass_table = mass_table

    def load_candidates(self, psm_list):
        pass

    def _get_alternative_tags(self, mass, max_n):
        masses = self.mass_table[
            (self.mass_table.mass < (mass+self.tolerance)) &
            (self.mass_table.mass > (mass-self.tolerance))
        ]
        alternatives = list(chain(*masses.aa.tolist()))

        # Parse the tags to peptidoform format
        alternatives = [
            self.parse_tag_to_unimod(alternative) for alternative in alternatives
        ]

        if len(alternatives) < max_n:
            return alternatives
        
        return random.sample(alternatives, max_n)
    
    def generate_peptide_variants(
            self,
            peptidoform,
            n=100,
            peptide_evidence=None,
            n_terminal=False,
            c_terminal=False,
            ambiguous=False
        ):
        peptide_variants = []

        if ambiguous and peptide_evidence is not None:
            peptide_variants += self._get_ambiguous_variants(
                peptidoform, peptide_evidence, n=n
            )

        if n_terminal:
            peptide_variants += self._get_nterm_variants(
                peptidoform, n=n
            )

        if c_terminal:
            peptide_variants += self._get_cterm_variants(
                peptidoform, n=n
            )

    def _get_nterm_tags(self, peptidoform, n, n_amino_acids=3):
        variants = []
        
        n_amino_acids = min(n_amino_acids, 5)

        terminal_mod = peptidoform.properties.get('n_term', None)
        tag = peptidoform[:n_amino_acids]
        indices = [0, n_amino_acids]
        
        _, mass_tag = self.parse_tag_to_pt(
            tag=tag,
            n_terminal_mod=terminal_mod
        )

        # Generate variants which have an n-terminal modification
        if  isinstance(terminal_mod, list):
            mass_mod = ptmass.calculate_mass(self.aa_comp[self.unimod2pyteomics[str(terminal_mod[0])]+"-"] - {'H': 1})
            mass_tag_no_mod = mass_tag - mass_mod
            alternatives_same_mod = self._get_alternative_tags(mass=mass_tag_no_mod, max_n=n)
            for alternative in alternatives_same_mod:
                alternative.properties['n_term'] = terminal_mod
            variants += alternatives_same_mod

        variants += self._get_alternative_tags(mass=mass_tag, max_n=n)

        if len(variants) > n:
            variants = random.sample(variants, n)
    
        return variants, indices

    def _get_cterm_tags(self, peptidoform, n, n_amino_acids=3):
        variants = []

        tag = peptidoform[-n_amino_acids:]
        indices = [-n_amino_acids, len(peptidoform)]
        _, mass_tag = self.parse_tag_to_pt(
            tag=tag
        )

        variants += self._get_alternative_tags(mass=mass_tag, max_n=n)

        if len(variants) > n:
            variants = random.sample(variants, n)
    
        return variants, indices

    def get_terminal_variants(self, peptidoform, n, n_amino_acids, terminus):
        if terminus == 'n':
            tag_variants, indices = self._get_nterm_tags(peptidoform, n, n_amino_acids)

        elif terminus == 'c':
            tag_variants, indices = self._get_cterm_tags(peptidoform, n, n_amino_acids)
    
        else:
            raise ValueError(f"Terminus can only be 'n' or 'c'. '{terminus}' was provided.")

        return [
            self.replace_tag(peptidoform=peptidoform, tag=tag, indices=indices)
            for tag in tag_variants
        ]

    def get_ambiguous_variants(self, peptidoform, peptide_evidence: PeptideEvidence, n):
        
        peptidoform_variants = []
        tag_indices = []
        ambiguous_sites = []
        
        # Generate tags for each ambiguous position
        for i, indices in enumerate(peptide_evidence.get_ambiguous_tag_idx()):
            tag = peptidoform[indices[0]: indices[1]]

            # Important to take into account potential n-term modifications
            if indices[0]==0:
                variants, _ = self._get_nterm_tags(
                    peptidoform=peptidoform,
                    n=n,
                    n_amino_acids=indices[1]
                )
            
            else:
                variants = []
                if indices[1] - indices[0] > 5:
                    idx0 = indices[0]
                    idx1 = indices[0]+5
                    for _ in range(indices[1]-idx1):
                        tag = peptidoform[idx0:idx1]
                        _, mass_tag = self.parse_tag_to_pt(
                            tag=tag
                        )
                        variants += self._get_alternative_tags(mass=mass_tag, max_n=n)
                        idx0 += 1
                        idx1 += 1

                else:
                    tag = peptidoform[indices[0]: indices[1]]
                    _, mass_tag = self.parse_tag_to_pt(
                        tag=tag
                    )

                    variants += self._get_alternative_tags(mass=mass_tag, max_n=n)

            ambiguous_sites.append(variants)
            tag_indices.append(indices)


        # Collect n combinations of the tag variants
        # all_combinations = list(product(*ambiguous_sites))
        # n = min(n, len(all_combinations))
        # selected_combinations = random.sample(all_combinations, n)
        selected_combinations = []
        for _ in range(n):
            combination = [random.choice(site) for site in ambiguous_sites]
            selected_combinations.append(combination)

        # Replace the tags at the positions matching the indices
        for combination in selected_combinations:
            modified_peptidoform = peptidoform
            current_indices = tag_indices.copy()
            offset = 0  # Tracks the index shift caused by replacements

            for tag, indices in zip(combination, current_indices):
                start, end = indices
                start += offset
                end += offset

                # Replace the tag and calculate the offset change
                modified_peptidoform = self.replace_tag(
                    peptidoform=modified_peptidoform,
                    tag=tag,
                    indices=[start, end]
                )
                offset += len(tag.parsed_sequence) - (end - start)

            peptidoform_variants.append(modified_peptidoform)
        
        peptidoform_variants = list(set(peptidoform_variants))
        return peptidoform_variants


    def parse_tag_to_pt(self, tag, n_terminal_mod=None):
        """
        Parse a tag to pyteomics compatible format.
        
        Parameters
        ----------
        tag: list
            List of two-sized tuples (amino acid, modification).
        n_terminal_mod: list or None
            The N-terminal modification.

        Returns
        -------
        tuple
            str: Reconstructed tag readible with psm_utils
            float: The mass of the tag
        """
        tag_peptide_reconstr = ""
        tag_mass = ""
        for aa, mods in tag:
            tag_peptide_reconstr+=aa
            if mods:
                for mod in mods:
                    mod = str(mod)
                    tag_peptide_reconstr+="["+mod+"]"
                    tag_mass+=self.unimod2pyteomics[mod]
            tag_mass += aa

        if isinstance(n_terminal_mod, list):
            n_term_mod_str = str(n_terminal_mod[0])
            tag_peptide_reconstr = f"[{n_term_mod_str}]-" + tag_peptide_reconstr
            tag_mass = self.unimod2pyteomics[n_term_mod_str] + "-" + tag_mass + "-OH"
        
        return (tag_peptide_reconstr, ptmass.calculate_mass(tag_mass, aa_comp=self.aa_comp))    

    def parse_tag_to_unimod(self, tag):
        for pt, unimod in self.pyteomics2unimod.items():
            tag = tag.replace(pt, unimod)
        return Peptidoform(tag)
    
    def replace_tag(self, peptidoform, tag, indices):
        start, end = indices
        properties, parsed_sequence = peptidoform.properties.copy(), [i for i in peptidoform.parsed_sequence]
        parsed_sequence[start: end] = tag.parsed_sequence
        if start == 0:
            properties['n_term'] = tag.properties['n_term']
    
        return Peptidoform(proforma.to_proforma(sequence=parsed_sequence, **properties))


UNIMOD2PYTEOMICS = {
    "UNIMOD:4": "c",
    "UNIMOD:35": "ox",
    "UNIMOD:1": "Ac",
    "UNIMOD:5": "Carbamyl",
    "UNIMOD:385": "Ammonialoss",
    "+25.980265": "CarbamylAmmonialoss",
    "UNIMOD:7": "deam"
}
PYTEOMICS2UNIMOD = {
    "cC": "C[UNIMOD:4]",
    "oxM": "M[UNIMOD:35]",
    "Ac-": "[UNIMOD:1]-",
    "Carbamyl-": "[UNIMOD:5]-",
    "Ammonialoss-": "[UNIMOD:385]-",
    "CarbamylAmmonialoss-": '[+25.980265]-',
    "deamQ": 'Q[UNIMOD:7]',
    "deamN": 'N[UNIMOD:7]'
}

def get_aa_composition():
    db = ptmass.Unimod()
    aa_comp = ptmass.std_aa_comp
    aa_comp["ox"] = db.by_id(35)["composition"]
    aa_comp["c"] = db.by_id(4)["composition"]
    aa_comp['ac'] = db.by_title('Acetyl')['composition']
    aa_comp['Ac-'] = aa_comp['ac'] + {'H': 1}
    aa_comp['deam'] = db.by_id(7)['composition']
    aa_comp['Ammonialoss-'] = db.by_id(385)["composition"] + {'H':1}
    aa_comp['Carbamyl-'] = db.by_id(5)["composition"] + {'H':1}
    aa_comp['CarbamylAmmonialoss-'] = db.by_id(5)["composition"] + db.by_id(385)["composition"] + {'H': 1}
    return aa_comp

def generate_mass_table(aa_comp=get_aa_composition()):
    # Create the mass dataframe
    aas = [i for i in "ADEFGHLKMNPQRSTVWY"] + ["cC", "oxM"]
    aa_2_comp = list(["".join(x) for x in product(aas, repeat=2)])
    aa_3_comp = list(["".join(x) for x in product(aas, repeat=3)])
    aa_4_comp = list(["".join(x) for x in product(aas, repeat=4)])
    aa_5_comp = list(["".join(x) for x in product(aas, repeat=5)])

    aa_comps = aa_2_comp+aa_3_comp+aa_4_comp+aa_5_comp

    aa_df = pd.DataFrame(aa_comps).rename(columns={0: "aa"})
    aa_df["mass"] = aa_df.progress_apply(lambda x: np.round(
            ptmass.calculate_mass(x["aa"], aa_comp=aa_comp),
            decimals=6
        ), axis=1)
    mass_to_aa = aa_df.groupby("mass").aa.apply(list).reset_index().sort_values("mass")
    return mass_to_aa


def candidatelist_to_psmlist(candidates, psm):
    psm_list = []
    for candidate in candidates:
        new_psm = psm.model_copy(deep=True)
        new_psm.peptidoform = candidate
        psm_list.append(new_psm)
    
    return PSMList(psm_list=psm_list)

def generate_candidate_set(
        spectrum,
        cgenerator: CandidateGenerator
    ):

    psm_gt = PSM(
        peptidoform=spectrum.psm_gt.peptidoform,
        spectrum_id=spectrum.spectrum_id,
        **spectrum.properties['properties']
    )

    peptide_evidences = [spectrum.psm_gt.peptide_evidence] 
    for psm in spectrum.psm_candidates:
        if psm.evaluation['score_ms2rescore_features'].error_type == 'match':
            continue
        else:
            peptide_evidences.append(psm.peptide_evidence)
    engine_peptides = [pe.peptidoform for pe in peptide_evidences]


    ambiguous = []
    terminal_changes = []
    # GENERATE PEPTIDE VARIANTS
    for peptide_evidence in peptide_evidences:

        # most difficult cases ?
        ambiguous += cgenerator.get_ambiguous_variants(
            peptidoform=peptide_evidence.peptidoform,
            peptide_evidence=peptide_evidence,
            n=100
        )

        # Less difficult
        terminal_changes += cgenerator.get_terminal_variants(
            peptide_evidence.peptidoform,
            n=100,
            n_amino_acids=3,
            terminus='n'
        )
        terminal_changes += cgenerator.get_terminal_variants(
            peptide_evidence.peptidoform,
            n=100,
            n_amino_acids=3,
            terminus='c'
        )

        # # Crappy cases
        # terminal_changes += cgenerator.get_terminal_variants(
        #     peptide_evidence.peptidoform,
        #     n=500,
        #     n_amino_acids=4,
        #     terminus='n'
        # )
        # terminal_changes += cgenerator.get_terminal_variants(
        #     peptide_evidence.peptidoform,
        #     n=500,
        #     n_amino_acids=4,
        #     terminus='c'
        # )

    # Ensure no duplicates
    candidates = list(set(terminal_changes))
    # Ensure ground-truth is not represented in candidate list
    candidates = [candidate for candidate in candidates if candidate not in engine_peptides]

    candidate_psmlist = candidatelist_to_psmlist(candidates, psm_gt)
    candidate_psmlist_ambiguous = candidatelist_to_psmlist(ambiguous, psm_gt)

    # Annotate the variants
    for psm in candidate_psmlist:
        psm['source'] = 'Terminal-variant'
    for psm in candidate_psmlist_ambiguous:
        psm['source'] = "Ambiguous"

    # Limit to 1000 PSMs.
    if len(candidate_psmlist) > 1000:

        candidate_psmlist = PSMList(psm_list=random.sample(list(candidate_psmlist), 1000))
    candidate_psmlist = candidate_psmlist + candidate_psmlist_ambiguous
    candidate_psmlist['is_decoy'] = [True]*(len(candidate_psmlist)-1) + [False]
    return candidate_psmlist