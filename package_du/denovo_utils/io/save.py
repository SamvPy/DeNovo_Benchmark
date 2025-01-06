import pickle
import os
from psm_utils import PSMList
from psm_utils.io import write_file

def save_pickle(data, path):
    with open(path, 'wb') as f:
        pickle.dump(data, f)

def save_mokapot_models(models, save_folder):
    if save_folder is not None:
        for i, model in enumerate(models):
            model.save(
                os.path.join(save_folder, f"mokapot_model_{i}.pkl")
            )

def save_features(psm_list, save_path):
    psm_features = psm_list.to_dataframe()[["spectrum_id", "rank", "run", "source", "rescoring_features"]]
    psm_features.to_parquet(save_path)

def save_psmlist(psm_list, save_path):
    '''
    After saving psm_list, the original PSMList object has his rescoring features removed.
    '''
    if isinstance(psm_list, PSMList):
        # Fix storage of a peptide_evidence object.
        if "peptide_evidence" in psm_list[0]['metadata'].keys():
            for psm in psm_list:
                # Values in the metadata dict must be stringified according to parquet schema
                # When reading, these must be converted again
                metadata = psm['metadata']
                evidence = str(list(metadata['peptide_evidence'].evidence))
                evidence_labels = str(list(metadata['peptide_evidence'].evidence_labels))
                metadata['PE_evidence'] = evidence
                metadata['PE_evidence_labels'] = evidence_labels
        
        for psm in psm_list:
            psm['metadata'] = {k: str(v) for k, v in psm['metadata'].items() if k!='peptide_evidence'}
            psm['provenance_data'] = {k: str(v) for k, v in psm['provenance_data'].items()}
            psm['rescoring_features'] = {}

        write_file(psm_list=psm_list, filename=save_path, filetype='parquet')
    else:
        raise Exception(f"Cannot write object with type {type(psm_list)} to parquet. Must be type {type(PSMList)}.")

