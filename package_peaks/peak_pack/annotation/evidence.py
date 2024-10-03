from psm_utils import Peptidoform

class PeptideEvidence:
    def __init__(self, peptidoform, ion_matrix, evidence_labels):

        if isinstance(peptidoform, str):
            self.peptidoform = Peptidoform(peptidoform)
        else:
            self.peptidoform = peptidoform

        self.evidence_labels = evidence_labels
        self.evidence = ion_matrix.any(axis=0)
        self.ambiguous_tag_indices = self.get_ambiguous_tags()    

    def __repr__(self):

        peptide_repr = ""
        first_tags = [i[0] for i in self.ambiguous_tag_indices]
        end_tags = [i[1]-1 for i in self.ambiguous_tag_indices]

        for i, aa in enumerate(self.peptidoform.parsed_sequence):
            mod = ""
            if aa[1]:
                mod = "[" + str(aa[1][0]) + "]"
            aa_repr = aa[0]+mod

            if i in first_tags:
                peptide_repr += "<"
            peptide_repr += aa_repr
            if i in end_tags:
                peptide_repr += ">"
        peptide_repr += "/"+str(self.peptidoform.precursor_charge)
        return peptide_repr
    
    def get_ambiguous_tags(self):
        """
        Creates a list of tuples (start_index, end_index) of amino acids which can be substituted with anything with equal mass
        """
        prev_evidence = True
        isobaric_parts = []
        isobaric_part = []
        for i, aa in enumerate(self.peptidoform.parsed_sequence):

            if i == len(self.peptidoform.parsed_sequence)-1:
                break
            
            evidence = self.evidence[i]
            
            if prev_evidence:
                if not evidence:
                    isobaric_part.append(i)
            if not prev_evidence:
                if evidence:
                    isobaric_part.append(i+1)
                    isobaric_parts.append(isobaric_part)
                    isobaric_part = []
            prev_evidence=evidence

        if not prev_evidence:
            isobaric_part.append(i+1)
            isobaric_parts.append(isobaric_part)
        return isobaric_parts