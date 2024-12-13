import pickle
import pandas as pd
import logging
from glob import glob
import json
import os
from psm_utils.io import read_file
from psm_utils import PSMList
from ms2rescore.feature_generators import FEATURE_GENERATORS
from ms2rescore.feature_generators.deeplc import DeepLCFeatureGenerator
from ms2rescore import parse_configurations

# Handle this somewhere else ? Risk of circular imports
from peak_pack.annotation.evidence import PeptideEvidence


logger = logging.getLogger(__name__)

def read_psmlist(path):
    with open(path, 'rb') as f:
        psmlist = pickle.load(f)
    return psmlist
    
def read_features(path):
    return pd.read_parquet(path)

def load_pickle(path):
    if not os.path.exists(path=path):
        raise Exception(f'Tried unpickling but file not found: {path}')

    with open(path, 'rb') as f:
        data = pickle.load(f)
    return data

def load_deeplc(
        model_path_folder: str,
        calibration_path: str
    ) -> DeepLCFeatureGenerator:
    logging.info("Loading deeplc from cache")
    model_paths = []
    for p in glob(os.path.join(model_path_folder, "*.keras")):
        model_paths.append(p)
    
    calibration_params = load_pickle(path=calibration_path)
    deeplc_kwargs = {"path_model": model_paths}
    deeplc_kwargs.update(calibration_params)
    deeplc_fgen = FEATURE_GENERATORS["deeplc"](**deeplc_kwargs)
    deeplc_fgen.deeplc_predictor.calibrate_dict = deeplc_kwargs['calibrate_dict']
    deeplc_fgen.deeplc_predictor.calibrate_min = deeplc_kwargs['calibrate_min']
    deeplc_fgen.deeplc_predictor.calibrate_max = deeplc_kwargs['calibrate_max']

    logging.info(f"Loaded DeepLC models: {deeplc_fgen.deeplc_predictor.model}")
    logging.info(f"Loaded calibration sets at {calibration_path}")
    return deeplc_fgen

def load_configuration(
        path_config_ms2rescore,
        psm_file,
        spectrum_path
    ):
    with open(path_config_ms2rescore) as config_file:
        configuration = json.load(config_file)
    configuration["ms2rescore"]["psm_file"] = psm_file
    configuration["ms2rescore"]["spectrum_path"] = spectrum_path
    configuration = parse_configurations(configurations=configuration)
    return configuration

def load_mokapot(mokapot_folder):
    models = []
    for i in range(3):
        with open(os.path.join(mokapot_folder, f"mokapot_model_{i}.pkl"), "rb") as f:
            models.append(pickle.load(f))
    logger.info(f"Loaded mokapot models from {mokapot_folder}")
    return models

def load_features(feature_path) -> dict:
    """
    Load a feature file as a dictionary where:
    - keys: spectrum_id
    - values: feature dictionary
    """
    features = pd.read_parquet(feature_path)
    if "rescoring_features" in features.columns:
        features = pd.concat(
            [
                features[features.columns[:3]].reset_index(drop=True),
                pd.DataFrame(features['rescoring_features'].tolist())
            ],
            axis=1
        )
    recs = {}
    # Iterate through the rows of the DataFrame
    for _, row in features.iterrows():
        spec_id = row["spectrum_id"]
        rank = row["rank"]
        # Convert feature columns to a dictionary for the current row
        feature_dict = row[features.columns[4:]].to_dict()
        
        # Ensure there is a dictionary for this spec_id
        if spec_id not in recs:
            recs[spec_id] = {}
        
        # Store the feature_dict under the corresponding rank
        recs[spec_id][rank] = feature_dict
    return recs

def load_psmlist(psm_path) -> PSMList:
    psm_list = read_file(
        psm_path,
        filetype='parquet'
    )

    # Reload the peptide_evidence objects
    if ('PE_evidence_labels' in psm_list[0]['metadata'].keys() and
        'PE_evidence' in psm_list[0]['metadata'].keys()):
        for psm in psm_list:
            metadata = {k: eval(v) for k, v in psm['metadata'].items()}
            metadata['peptide_evidence'] = PeptideEvidence(
                peptidoform=psm['peptidoform'],
                evidence_labels=eval(psm['metadata']['PE_evidence_labels']),
                evidence=eval(psm['metadata']['PE_evidence'])
            )
            provenance_data = {k: eval(v) for k, v in psm['provenance_data'].items()}

            psm['metadata'] = metadata
            psm['provenance_data'] = provenance_data

    return psm_list

def load_psmlist_and_features(psm_path, feature_path):
    psm_list = load_psmlist(psm_path)
    features = load_features(feature_path)

    psm_list['rescoring_features'] = [{} for i, _ in enumerate(psm_list)]
    _ = [
        psm["rescoring_features"].update(
            features[
                psm["spectrum_id"] # Grouped by spectrum_id
            ][
                psm['rank'] # Multiple ranks possible per spectrum_id
            ]
        ) 
        for psm 
        in psm_list
    ]

    return psm_list
