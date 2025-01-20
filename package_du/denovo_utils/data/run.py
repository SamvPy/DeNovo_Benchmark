from typing import Optional, List
from peak_pack.annotation import PeptideEvidence
from psm_utils import Peptidoform
from tqdm import tqdm
from .psm import PSM
from .spectrum import Spectrum
from ..parsers import DenovoEngineConverter
import logging
import pandas as pd

logger = logging.getLogger(__name__)

ENGINE_NAME_MAPPER = {
    "Casanovo4.2.0": "casanovo",
    "ContraNovo": "contranovo",
    "InstaNovo": "instanovo",
    "NovoB": "novob",
    "PepNet": "pepnet"
}

class Run:
    def __init__(self, run_id, ignore_IL=True):
        self.run_id = run_id
        self.spectra: dict[str, Spectrum] = {}  # Dictionary to map spectrum_id to Spectrum object
        self.ignore_IL = ignore_IL
    
    def __repr__(self):
        engines = []
        for spectra in self.spectra.values():
            for psm in spectra.psm_candidates:
                if psm.engine_name not in engines:
                    engines.append(psm.engine_name)
        return f"{len(self.spectra)} spectra loaded from {len(engines)} engines ({engines})."

    def load_data(self, psmlist, score_names, is_ground_truth=False, filter_by_gt=True):

        for psm_ in tqdm(psmlist):
            # Only ground truth psms can add a Spectrum object to the run.
            if is_ground_truth:
                if psm_['is_decoy'] or psm_['qvalue']>.01:
                    continue
                properties = {
                    "precursor_mz": psm_["precursor_mz"],
                    "retention_time": psm_['retention_time'],
                    "ion_mobility": psm_['ion_mobility'],
                    "is_decoy": psm_['is_decoy'],
                }
                spectrum = Spectrum(psm_['spectrum_id'], properties=properties)
                self.add_spectrum(spectrum=spectrum)

            spectrum = self.get_spectrum(psm_["spectrum_id"])
            
            if not isinstance(spectrum, Spectrum):
                if filter_by_gt:
                    continue
                # Add spectra if psms without a gt still needs to be added.
                else:
                    spectrum = Spectrum(psm_["spectrum_id"])
                    self.add_spectrum(spectrum)
            
            engine_name = psm_["source"]
            if psm_["source"] in ENGINE_NAME_MAPPER.keys():
                engine_name = ENGINE_NAME_MAPPER[psm_["source"]]
            elif psm_['source'] is None:
                if is_ground_truth:
                    engine_name = "Ground-Truth"
                else:
                    engine_name = "Unspecified"
            else:
                engine_name = psm_['source']

            if self.ignore_IL:
                peptidoform = leucine_to_isoleucine(psm_["peptidoform"])
            else:
                peptidoform = psm_["peptidoform"]

            psm = PSM(
                peptidoform=peptidoform,
                score=psm_['score'],
                engine_name=engine_name,
                rank=psm_['rank'],
                aa_score=extract_metadata(psm_, 'aa_scores', infer_datatype=True),
                peptide_evidence=PeptideEvidence.load(
                    extract_metadata(psm_, 'peptide_evidence')
                )
            )
            
            # Add ms2rescore score if rescored.
            for score_name in score_names:
                if score_name not in psm_['provenance_data'].keys():
                    continue

                psm.scores.add_score(
                    score=psm_['provenance_data'][score_name],
                    metadata=score_name,
                    score_type='peptide'
                )

            spectrum.add_psm(psm, is_ground_truth=is_ground_truth)


    def load_refinement(self, psmlist):
        for psm_ in psmlist.get_rank1_psms():
            spectrum = self.get_spectrum(psm_["spectrum_id"])

            # Cannot load refinements if spectrum has no associated PSMs
            if not isinstance(spectrum, Spectrum):
                continue


            if self.ignore_IL:
                peptidoform = leucine_to_isoleucine(psm_["peptidoform"])
            else:
                peptidoform = psm_["peptidoform"]

            # Parse refinement into a PSM
            psm = PSM(
                peptidoform=peptidoform,
                score=psm_['score'],
                engine_name=psm_['source'],
                rank=psm_['rank'],
                aa_score=extract_metadata(psm_, 'aa_scores', infer_datatype=True),
                peptide_evidence=PeptideEvidence.load(
                    extract_metadata(psm_, 'peptide_evidence')
                )
            )
            # Extract base PSM to load the refinement into
            base_prediction = psm_["metadata"]["base_prediction"]
            if base_prediction in ENGINE_NAME_MAPPER.keys():
                base_prediction = ENGINE_NAME_MAPPER[base_prediction]

            psm_base_list = spectrum.get_psms_by_engine(base_prediction)
            if len(psm_base_list)==0:
                logging.warning(f"No base prediction ({base_prediction}-{psm_['source']}) for spectrum: {psm_['spectrum_id']} ")
                continue
            elif len(psm_base_list) > 1:
                raise Exception(f"Cannot decide on base psm. {len(psm_base_list)} psms found: {psm_base_list}")
            else:
                psm_base = psm_base_list[0]

            # Add additional scoring if present
            if 'score_ms2rescore' in psm_['provenance_data'].keys():
                psm.scores.add_score(
                    score=psm_['provenance_data']['score_ms2rescore'],
                    metadata='score_ms2rescore',
                    score_type='peptide'
                )
        
            psm_base.add_refinement(psm)
            
            # Spectralis also rescores original identification. Add this to score object.
            if psm_["source"] == "Spectralis":
                try:
                    init_score = eval(psm_['provenance_data']['init_score_spectralis'])
                except TypeError:
                    init_score = psm_["provenance_data"]["init_score_spectralis"]
                psm_base.scores.add_score(
                    score=init_score,
                    metadata="Spectralis",
                    score_type="peptide",
                    overwrite=True
                )

    def add_spectrum(self, spectrum: Spectrum):
        if spectrum.spectrum_id in self.spectra.keys():
            raise Exception(f"Cannot add an existing spectrum. {spectrum.spectrum_id}")
        self.spectra[spectrum.spectrum_id] = spectrum
    
    def get_spectrum(self, spectrum_id) -> Spectrum:
        return self.spectra.get(spectrum_id, None)
    
    def get_correct_spectra(self, engines=[], correctness='match', score_metadata="ms2rescore", exclusive=False) -> 'Run':
        correctness_order = ['match', 'isobaric_aa', 'isobaric_peak', 'Higher score', "Lower score", "Different"]
        correctness_index = correctness_order.index(correctness)

        spectra_dict = {}

        for spectrum_id, spectra in self.spectra.items():
            psms_engine_specific = []
            psms_true = []
            for psm in spectra.psm_candidates:
                
                if correctness_order.index(psm.evaluation[score_metadata].error_type) <= correctness_index:
                    psms_true.append(psm)
                    
                    if psm.engine_name in engines:
                        psms_engine_specific.append(psm)
                        continue
            
            if len(psms_engine_specific) > 0:
                if exclusive and len(psms_true) > 1:
                    continue

                new_spectrum = Spectrum(spectrum_id)
                new_spectrum.add_psm(spectra.psm_gt, is_ground_truth=True)
                for psm in psms_engine_specific:
                    new_spectrum.add_psm(psm)
                spectra_dict[spectrum_id] = new_spectrum
            
        new_run = Run(self.run_id)
        new_run.spectra = spectra_dict
        return new_run

    def get_unique_spectra(self, engine) -> 'Run':
        spectra_dict = {}
        for spectrum_id, spectra in self.spectra.items():
            psms = []
            
            
            psm = spectra.get_psms_by_engine(engine_name=engine)
            if len(psm)==0:
                continue

            if len(psm)==len(spectra):
                spectra_dict[spectrum_id] = spectra
            
        new_run = Run(self.run_id)
        new_run.spectra = spectra
        return new_run

    def get_common_spectra(self, engines, include_gt=True) -> 'Run':
        spectra = {}
        for spectrum_id, spectrum in self.spectra.items():
            psms = []
            common = True
            for engine in engines:
                psm = spectrum.get_psms_by_engine(engine_name=engine)
                if len(psm) == 0:
                    common = False
                    break
                psms += psm

            if common:
                if include_gt:
                    psms += [spectrum.psm_gt]
                spectra[spectrum_id] = spectrum
        
        new_run = Run(self.run_id)
        new_run.spectra = spectra
        return new_run

    @property
    def modifications(self):
        mod_list = []
        for spectrum in self.spectra.values():
            mod_list += spectrum.psm_gt.modifications
            mod_list = list(set(mod_list))
        return [mod.name for mod in mod_list]

    def get_modified_spectra(self, modification, exclude_modifications=[]):
        mod_specid = {"": []}

        for spectrum in self.spectra.values():
            spec_mods = spectrum.psm_gt.modifications
            has_excl_mod=False
            for ex_mod in exclude_modifications:
                if ex_mod in [spec_mod.name for spec_mod in spec_mods]:
                    has_excl_mod=True
                    break
            if has_excl_mod:
                continue
    
            for mod in spec_mods:
                if mod.name not in mod_specid.keys():
                    mod_specid[mod.name] = [spectrum.spectrum_id]
                else:
                    mod_specid[mod.name].append(spectrum.spectrum_id)

            if spec_mods == []:
                mod_specid[""].append(spectrum.spectrum_id)

        if modification not in mod_specid.keys():
            print(f"No GT spectra with {modification}. Found modifications: {list(mod_specid.keys())}")
            return []
        
        return mod_specid[modification]
    
    def get_features(self, features_path: str, keep_spectrum_ids=False):

        spectrum_ids = list(self.spectra.keys())
        features = pd.read_parquet(features_path)
        features = features[features.spectrum_id.isin(spectrum_ids)]
        spectrum_ids_col = features['spectrum_id'].reset_index(drop=True)

        features = pd.DataFrame(features.rescoring_features.tolist())

        if keep_spectrum_ids:
            features["spectrum_id"] = spectrum_ids_col
        return features

def extract_metadata(psm, key, infer_datatype=False):
    if key in psm['metadata'].keys():
        el = psm['metadata'][key]
        if infer_datatype:
            try:
                return eval(el)
            except:
                return el
        return el
    return None

def leucine_to_isoleucine(peptidoform):
    return Peptidoform(peptidoform.proforma.replace("L", "I"))