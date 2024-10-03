import numpy as np
from ms2rescore.feature_generators.base import FeatureGeneratorBase
from tqdm import tqdm
from psm_utils import Peptidoform, PSMList
from ms2pip.core import annotate_spectra
from ..utils.annotation import AnnotatedSpectrum
from spectrum_utils import spectrum as sus
from itertools import chain
from typing import Optional
from pyteomics import mgf

def get_annotated_spectrum(
        mgf_file,
        psm,
        ions="abcxyzIp",
        neutral_losses=False
    ):
    mgf_spectrum = mgf_file.get_by_id(psm["spectrum_id"])
    spectrum = sus.MsmsSpectrum(
        identifier=mgf_spectrum["params"]["title"],
        precursor_mz=mgf_spectrum["params"]["pepmass"][0],
        precursor_charge=int(mgf_spectrum["params"]["charge"][0]),
        mz=mgf_spectrum["m/z array"],
        intensity=mgf_spectrum["intensity array"],
        retention_time=mgf_spectrum["params"]["rtinseconds"]
    ).annotate_proforma(psm["peptidoform"].proforma, 50, "ppm", ions, neutral_losses=neutral_losses, max_ion_charge=2)
    return spectrum

def get_missing_fragmentations(
        ion_vectors,
        evidence_required = ["b1", "b2", "y1", "y2"],
    ):
    """
    B: starts from N-terminus
    Y: starts from C-terminus
    """
    length_peptide = len(ion_vectors[evidence_required[0]])
    fragmentation_sites = {i: [] for i in range(length_peptide)} # Starts from N-term
    for evidence_type in evidence_required:

        
        for i, ion_intensity in enumerate(ion_vectors[evidence_type]):
            detected = ion_intensity>0
            if evidence_type.startswith("b"):
                fragmentation_sites[i].append(detected)

            elif evidence_type.startswith("y"):
                fragmentation_sites[length_peptide-i-1].append(detected)

            else:
                raise NotImplementedError(f"{evidence_type} not implemented.")
        
    
    fragmentation_site_detected = []
    for fragmentation_site in fragmentation_sites.values():
        fragmentation_site_detected.append(True in fragmentation_site)
    
    return np.array(fragmentation_site_detected)

def get_ambiguous_tags(peptidoform, evidence_bool):
    """
    Creates a list of tuples (start_index, end_index) of amino acids which can be substituted with anything with equal mass
    """
    prev_evidence = True
    isobaric_parts = []
    isobaric_part = []
    for i, aa in enumerate(peptidoform.parsed_sequence):

        if i == len(peptidoform.parsed_sequence)-1:
            break
        
        evidence = evidence_bool[i]
        
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


class PeptideEvidence:
    def __init__(self, proforma, evidence_bool):
        self.peptidoform = Peptidoform(proforma)
        self.evidence = evidence_bool
        self.ambiguous_tag_indices = get_ambiguous_tags(self.peptidoform, self.evidence)    

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
    
class MissingFragmentationSiteFGen(FeatureGeneratorBase):
    def __init__(
            self,
            *args,
            spectrum_path: Optional[str] = None,
            only_by_ions = True,
            **kwargs,
    ) -> None:
        super().__init__(*args, **kwargs)
        self.spectrum_path = spectrum_path
        self.mgf_file = mgf.read(self.spectrum_path)
        self.only_by_ions = only_by_ions
    
    @property
    def feature_names(self):
        [
            "missing_frag_count",
            "missing_frag_longest_sequence",
            "missing_frag_from_n",
            "missing_frag_from_c",
            "missing_frag_proportion"
        ]
    
    def add_features(self, psm_list: PSMList) -> None:
        psm_dict = psm_list.get_psm_dict()

        current_run = 1
        total_runs = sum(len(runs) for runs in psm_dict.values())
        for runs in psm_dict.values():
            for run, psms in runs.items():

                psm_list_run = PSMList(psm_list=list(chain.from_iterable(psms.values())))

                if self.only_by_ions:
                    pr_result = annotate_spectra(
                        psms=psm_list_run,
                        spectrum_file=self.spectrum_path,
                        model="HCDch2"
                    )
                    for psm, processing_result in zip(psm_list_run, pr_result):
                        pe, features = self._calculate_features(psm=psm, processing_result=processing_result)
                        psm["rescoring_features"].update(
                            features
                        )
                        psm["metadata"].update(
                            {"peptide_evidence": pe}
                        )

                else:
                    for i, psm in tqdm(enumerate(psm_list_run), total=len(psm_list_run)):

                        meta, features = self._calculate_features(psm=psm)
                        psm["rescoring_features"].update(
                            features
                        )
                        psm["metadata"].update(
                            meta
                        )


    def processing_result_to_peptide_evidence(self, processing_result):
        
        ion_vectors = {}
        null_value = np.float32(np.log2(.001))
        ion_type_mapping = {
            "b": "b1",
            "y": "y1",
            "b2": "b2",
            "y2": "y2"
        }
        for ion_type, intensities in processing_result.observed_intensity.items():
            intensity_vector = []
            for intensity in intensities:
                if intensity == null_value:
                    intensity_vector.append(0.)
                else:
                    intensity_vector.append(1.)
            ion_vectors[ion_type_mapping[ion_type]] = intensity_vector

        return PeptideEvidence(
            processing_result.psm.peptidoform.proforma,
            evidence_bool=get_missing_fragmentations(ion_vectors=ion_vectors)
        )

    def peptide_evidence_to_features(self, psm, peptide_evidence, set_type=""):
        if peptide_evidence.ambiguous_tag_indices == []:
            return {
                f"{set_type}missing_frag_count": 0,
                f"{set_type}missing_frag_longest_sequence": 0,
                f"{set_type}missing_frag_from_n": -1,
                f"{set_type}missing_frag_from_c": -1,
                f"{set_type}missing_frag_proportion": 0
            }

        lengths_missing_frags = [
            i[1]-i[0]-1 for i in peptide_evidence.ambiguous_tag_indices
        ]
        missing_frag_count = len(peptide_evidence.ambiguous_tag_indices)
        missing_frag_longest_sequence = max(lengths_missing_frags)
        missing_frag_from_n = peptide_evidence.ambiguous_tag_indices[0][0]
        missing_frag_from_c = len(psm["peptidoform"])-peptide_evidence.ambiguous_tag_indices[-1][1]
        missing_frag_proportion = sum(lengths_missing_frags) / (len(psm["peptidoform"])-1)
        return {
            f"{set_type}missing_frag_count": missing_frag_count,
            f"{set_type}missing_frag_longest_sequence": missing_frag_longest_sequence,
            f"{set_type}missing_frag_from_n": missing_frag_from_n,
            f"{set_type}missing_frag_from_c": missing_frag_from_c,
            f"{set_type}missing_frag_proportion": missing_frag_proportion
        }

    def _calculate_features(self, psm, processing_result=None, ions="by"):

        if self.only_by_ions:
            peptide_evidence = self.processing_result_to_peptide_evidence(processing_result)
            features = self.peptide_evidence_to_features(psm, peptide_evidence=peptide_evidence)

        else:
            sa = AnnotatedSpectrum(
                spectrum=get_annotated_spectrum(mgf_file=self.mgf_file, psm=psm, ions=ions, neutral_losses=True),
                peplen=len(psm["peptidoform"].sequence),
                process=True
            )
            peptide_evidence = PeptideEvidence(
                sa.spectrum.proforma,
                evidence_bool=get_missing_fragmentations(sa.ion_vectors)
            )

            nl_bool = (
                get_missing_fragmentations(sa.nl_ion_vectors["-H2O"]) | 
                get_missing_fragmentations(sa.nl_ion_vectors["-NH3"])
            )
            peptide_evidence_nl = PeptideEvidence(
                sa.spectrum.proforma,
                evidence_bool=get_missing_fragmentations(sa.ion_vectors) | nl_bool
            )
            features = self.peptide_evidence_to_features(psm, peptide_evidence=peptide_evidence)
            features_nl = self.peptide_evidence_to_features(psm, peptide_evidence=peptide_evidence_nl, set_type="NL_")
            features.update(features_nl)
        
        metadata = {
            "annotated_spectrum": sa,
            "peptide_evidence": peptide_evidence
        }

        return metadata, features