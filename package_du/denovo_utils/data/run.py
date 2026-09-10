from typing import Optional, List
from peak_pack.annotation import PeptideEvidence
from psm_utils import Peptidoform, PSMList
from tqdm import tqdm
from .psm import PSM
from .spectrum import Spectrum
from ..parsers import DenovoEngineConverter
from ..parsers.constants import EXTENSIONS
from ..parsers.converters.spectralis import load_spectralis_rescoring
from ..io.read import load_psmlist
import logging
import pandas as pd
import os
import numpy as np

logger = logging.getLogger(__name__)

ENGINE_NAME_MAPPER = {
    "Casanovo4.2.0": "casanovo",
    "ContraNovo": "contranovo",
    "InstaNovo": "instanovo",
    "NovoB": "novob",
    "PepNet": "pepnet",
    "pihelixnovo": "pi-HelixNovo",
    "piprimenovo": "pi-PrimeNovo",
    "adanovo": "AdaNovo"
}

ENGINE_NAME_MAPPER_REVERSE = {
    'casanovo': 'Casanovo4.2.0',
    'contranovo': 'ContraNovo',
    'instanovo': 'InstaNovo',
    'novob': 'NovoB',
    'pepnet': 'PepNet',
    "pihelixnovo": "pi-HelixNovo",
    "piprimenovo": "pi-PrimeNovo",
    "adanovo": "AdaNovo"
}

tqdm.pandas()


class Run:
    def __init__(self, run_id):
        self.run_id = run_id
        self.spectra: dict[str, Spectrum] = {}
        self.loaded_goldstandard = False
        self.loaded_engines = {}
    
    def __repr__(self):
        return f"{self.run_id} ({len(self.spectra)} spectra)"
    
    def __len__(self):
        return len(self.spectra)

    def load_gold_standard(
            self,
            psmlist: PSMList,
            score_name=None,
            metadata_label: Optional[str]=None,
            engine_name=None,
            filter_on_qvalue: Optional[float]=None,
            filter_decoys: bool=True
        ):
        # Filtering search results to load gold standard
        length_psmlist = len(psmlist)
        if filter_on_qvalue:
            psmlist = psmlist[psmlist['qvalue']<filter_on_qvalue]
            print('Q-value ({}) filtering: {} PSMs filtered.'.format(
                filter_on_qvalue,
                length_psmlist - len(psmlist)
            ))
            length_psmlist = len(psmlist)
        if filter_decoys:
            psmlist = psmlist[~psmlist['is_decoy']]
            print('Decoy filtering: {} PSMs filtered.'.format(
                length_psmlist - len(psmlist)
            ))
            length_psmlist = len(psmlist)

        # Iterate over the PSMs
        for psm_ in tqdm(psmlist, desc=f'Loading PSMs in Run object ({self.run_id})'):
            spectrum_id = psm_['spectrum_id']

            # Create new spectrum entry if spectrum not in run object yet
            spectrum = self.get_spectrum(spectrum_id)
            if not isinstance(spectrum, Spectrum):
                spectrum = Spectrum(
                    spectrum_id=spectrum_id,
                    properties={
                        "precursor_mz": psm_["precursor_mz"],
                        "retention_time": psm_['retention_time'],
                        "ion_mobility": psm_['ion_mobility'],
                    }
                )
                self.add_spectrum(spectrum)
            
            # Create new, more versatile PSM object

            psm = PSM(
                peptidoform=psm_['peptidoform'],
                score=psm_['score'],
                engine_name=psm_['source'],
                rank=psm_['rank'],
                score_name=score_name,
                aa_score=None,
                is_ground_truth=True,
                peptide_evidence=PeptideEvidence.load(
                    extract_metadata(
                        psm_,
                        'peptide_evidence'
                    )
                )
            )
            # Label the ground-truth with a special label
            if metadata_label:
                psm.metadata[metadata_label] = psm_['metadata'][metadata_label]
            self.spectra[spectrum_id].add_psm(psm=psm, is_gold_standard=True)
    
    def load_psmlist(
            self,
            psmlist,
            score_names=['score_ms2rescore'],
            only_gold_standard_spectra=True,
    ):
        if only_gold_standard_spectra:
            print('Filtering out PSMS from spectra not identified with gold standard.')
            selection_bool = np.isin(psmlist['spectrum_id'], np.array(list(self.spectra.keys())))
            psmlist = psmlist[selection_bool]
        
        for psm_ in tqdm(psmlist, desc=f'Loading PSMs in Run object ({self.run_id})'):

            spectrum = self.get_spectrum(psm_['spectrum_id'])
            if not isinstance(spectrum, Spectrum) and not only_gold_standard_spectra:
                spectrum = Spectrum(psm_["spectrum_id"])
                self.add_spectrum(spectrum)

            psm = PSM(
                peptidoform=psm_['peptidoform'],
                score=psm_['score'],
                engine_name=psm_['source'],
                rank=psm_['rank'],
                is_ground_truth=False,
                aa_score=extract_metadata(psm_, 'aa_scores', infer_datatype=True),
                peptide_evidence=PeptideEvidence.load(
                    extract_metadata(psm_, 'peptide_evidence')
                )
            )
            for score_name in score_names:
                if score_name not in psm_['provenance_data'].keys():
                    continue
                psm.scores.add_score(
                    score=psm_['provenance_data'][score_name],
                    metadata=score_name,
                    score_type='peptide'
                )
            spectrum.add_psm(psm, is_gold_standard=False)

    def load_psmlist_refinement(
            self,
            psmlist,
            overwrite=False
    ):
        for psm_ in psmlist:
            spectrum = self.get_spectrum(psm_['spectrum_id'])
            if not isinstance(spectrum, Spectrum):
                continue

            psm = PSM(
                peptidoform=psm_['peptidoform'],
                score=psm_['score'],
                engine_name=psm_['source'],
                rank=psm_['rank'],
                is_ground_truth=False,
                aa_score=extract_metadata(psm_, 'aa_scores', infer_datatype=True),
                peptide_evidence=PeptideEvidence.load(
                    extract_metadata(psm_, 'peptide_evidence')
                )
            )
            # Extract base PSM to load refinement into
            base_prediction = psm_["metadata"]["base_prediction"]
            if base_prediction in ENGINE_NAME_MAPPER_REVERSE.keys():
                base_prediction = ENGINE_NAME_MAPPER_REVERSE[base_prediction]
            
            psm_base_list = spectrum.get_psms_by_engine(base_prediction)
            if len(psm_base_list)==0:
                logging.warning(f"No base prediction ({base_prediction}-{psm_['source']}) for spectrum: {psm_['spectrum_id']} ")
                continue
            elif len(psm_base_list) > 1:
                base_found = False
                for psm_base in psm_base_list:
                    if psm_base.rank==psm.rank:
                        base_found = True
                        break
                if not base_found:
                    logging.warning(f"Base PSM with rank {psm.rank} not loaded in object. spectrum_id: {spectrum.spectrum_id}")
                    continue
                    # raise Exception(f"Cannot decide on base psm. {len(psm_base_list)} psms found for spectrum: {spectrum.spectrum_id}")
            else:
                psm_base = psm_base_list[0]
            
            # Add additional scoring if present
            if 'score_ms2rescore' in psm_['provenance_data'].keys():
                psm.scores.add_score(
                    score=psm_['provenance_data']['score_ms2rescore'],
                    metadata='score_ms2rescore',
                    score_type='peptide',
                    overwrite=overwrite
                )
            psm_base.add_refinement(psm, overwrite=overwrite)

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
                    overwrite=overwrite
                )

    def get_spectrum(self, spectrum_id) -> Spectrum:
        return self.spectra.get(spectrum_id, None)
    
    def add_spectrum(self, spectrum: Spectrum):
        if spectrum.spectrum_id in self.spectra.keys():
            raise Exception(f"Cannot add an existing spectrum. {spectrum.spectrum_id}")
        self.spectra[spectrum.spectrum_id] = spectrum


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

def read_runs(
        root_denovo_output: str,
        root_ground_truth: str,
        root_refinement: str,
        root_mgf: str,
        filenames: list,
        engine_names: list[str],
        refinement_names: list[str],
        format_ground_truth='sage',
        ignore_IL=True,
        rescored_db=True,
        rescored_denovo=True,
        rescored_refinement=True
) -> dict[str, Run]:
    runs = {}
    # Each raw file and its results constitutes a distinct Run object
    for run_name in filenames:
        mgf_path = os.path.join(root_mgf, run_name + ".mgf")

        logging.info(f'Loading run: {run_name}')
        run = Run(run_id=run_name, ignore_IL=ignore_IL)

        # Formatting is a bit different in this case
        if rescored_db:
            path_gt = os.path.join(root_ground_truth, run_name, 'psmlist', 'ground_truth.parquet')
            psmlist_gt = load_psmlist(path_gt)
            score_names = ['score_ms2rescore']
        else:
            path_gt = os.path.join(root_ground_truth, run_name + EXTENSIONS[format_ground_truth])
            parser = DenovoEngineConverter.select(format_ground_truth)
            psmlist_gt = parser.parse(
                result_path=path_gt,
                mgf_path=mgf_path
            )
            score_names = []

        logging.info(f'Loaded ground-truth from {path_gt}')
        run.load_data(
            psmlist=psmlist_gt.get_rank1_psms(),
            score_names=score_names,
            is_ground_truth=True
        )

        # Iterate over all de novo engines
        for engine_name in engine_names:
            logging.info(f"Loading model results: {engine_name}")
            parser = DenovoEngineConverter.select(engine_name)

            if rescored_denovo:
                path_denovo_file = os.path.join(
                    root_denovo_output,
                    run_name,
                    'psmlist',
                    f'{engine_name}.parquet'
                )
                psmlist_denovo = load_psmlist(path_denovo_file)
            else:
                path_denovo_file = os.path.join(
                    root_denovo_output,
                    engine_name,
                    run_name + f'.{engine_name}.some_extension'
                )
                psmlist_denovo = parser.parse(
                    result_path=path_denovo_file,
                    mgf_path=mgf_path
                )

            logging.info(f'Loaded denovo results from {path_denovo_file}')
            run.load_data(
                psmlist=psmlist_denovo,
                score_names=score_names,
                is_ground_truth=False
            )

            # Load refinements
            for refinement_name in refinement_names:

                # Format and paths will be different when rescored.
                if rescored_refinement:
                    refinement_path = os.path.join(
                        root_refinement,
                        run_name,
                        'psmlist',
                        f'{engine_name}.{refinement_name}.parquet'
                    )
                    
                    if os.path.exists(refinement_path):
                        logging.info(f'Loading refinement: {refinement_name} {engine_name} from {refinement_path}.')
                        psmlist_refinement = load_psmlist(refinement_path)
                        run.load_refinement(psmlist_refinement)

                else:
                    refinement_path = os.path.join(
                        root_refinement,
                        refinement_name,
                        engine_name
                    )
                    logging.info(f'Loading refinement: {refinement_name} {engine_name} from {refinement_path}.')
                    
                    # Spectralis results are split up by rank
                    if refinement_name == 'spectralis':
                        df_spectralis = load_spectralis_rescoring(
                            root_path=refinement_path,
                            filename=run_name
                        )
                        run.load_spectralis_rescoring(
                            df=df_spectralis
                        )

                    if refinement_name == 'instanovoplus':
                        parser = DenovoEngineConverter.select('instanovoplus')
                        psmlist_instanovoplus = parser.parse(
                            result_path=os.path.join(refinement_path, run_name+'.csv'),
                            mgf_path=mgf_path
                        )
                        run.load_refinement(
                            psmlist=psmlist_instanovoplus
                        )
        runs[run_name] = run

    return runs