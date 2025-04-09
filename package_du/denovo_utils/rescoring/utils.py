import os
from glob import glob
import logging

import numpy as np
from psm_utils import PSMList

from ..parsers.constants import EXTENSIONS

logger = logging.getLogger(__name__)

def already_trained(
        save_psm_path,
        save_feature_path,
        save_mokapot_paths
    ):

    for p in [save_psm_path, save_feature_path] + save_mokapot_paths:
        if not os.path.exists(p):
            return False
    return True

def parse_config_paths(
        ground_truth_folder,
        mgf_folder,
        denovo_folder,
        postprocessing_folder,
        engines,
        post_processors=["instanovoplus", "spectralis"],
        *args,
        **kwargs
):
    mgf_psm_dict = {}

    raw_extension = 'mgf'
    if len(glob(os.path.join(mgf_folder, '*.mgf'))) == 0:
        raw_extension = 'mzML'

    for mgf_path in glob(os.path.join(mgf_folder, f"*.{raw_extension}")):
        filename = os.path.basename(mgf_path).split(".")[0]
        psm_file_path = os.path.join(ground_truth_folder, filename + '.*')
        
        # Parse psmlists with any extension.
        paths = []
        for psm_file_path_ in glob(psm_file_path):
            paths.append(psm_file_path_)
        if len(paths)==0:
            logger.warning(f'While parsing paths: No file found at {psm_file_path}')
            continue
        if len(paths) > 1:
            # Take path with lowest number of characters.
            # Sometimes a file is generated with feature_names...
            psm_file_path = paths[0]
            for p in paths[1:]:
                if len(p) < len(psm_file_path):
                    psm_file_path = p
            paths[0] = psm_file_path
            logger.warning(f'Multiple files matching {psm_file_path}. Taking {p}.')
        psm_file_path = paths[0]

        if not os.path.exists(mgf_path) or not os.path.exists(psm_file_path):
            continue

        psm_dict = {
            "psm_file": psm_file_path,
            "spectrum_path": mgf_path
        }

        denovo_dict = {}
        for engine in engines:
            for post_processor in post_processors:
                for result_path in glob(
                    os.path.join(postprocessing_folder, post_processor)+f"/{engine}/{filename}*"
                ):
                    engine_label = f"{post_processor}__{engine}"
                    denovo_dict[engine_label] = result_path
            else:
                denovo_path = os.path.join(
                    denovo_folder, engine, filename + f".{engine}" + EXTENSIONS[engine]
                )
                if os.path.exists(denovo_path):
                    denovo_dict[engine] = os.path.join(
                        denovo_folder, engine, filename + f".{engine}" + EXTENSIONS[engine]
                    )

        if denovo_dict == {}:
            continue
        psm_dict["denovo"] = denovo_dict
        mgf_psm_dict[filename] = psm_dict
    return mgf_psm_dict


def setup_paths(output_folder, run_name, filename):
    """
    Set up all the paths necessary for rescoring.

    These include the following:
    - feature_path: Path to a file containing the features
    - psm_path: Path to the psmlist
    - mokapot_model_folder: Folder where the mokapot models are stored
    - mokapot_model_path: speaks for itself
    - calibration_folder: Folder where calibration files for the feature generators
    are stored

    Parameters
    ----------
    output_folder: str
        The root folder where all rescoring results are stored. 
        This folder will contain subfolders named by each run. Within this,
        4 folders are created:
        - features: Feature file
        - psmlist: PSMList
        - mokapot: Trained mokapot models
        - calibration: Feature generation calibration files
    
    run_name: str
        The name of the mass spec run undergoing rescoring
    
    filename: str
        The name of this specific rescoring attempt. Typically is 'ground_truth' or 
        the name of a search engine.
    """
    save_feature_path = "{}/{}/{}/{}.parquet".format(
        output_folder,
        run_name,
        "features",
        filename
    )
    save_psm_path = "{}/{}/{}/{}.parquet".format(
        output_folder,
        run_name,
        "psmlist",
        filename
    )

    mokapot_model_folder = "{}/{}/mokapot".format(
        output_folder,
        run_name
    )
    save_mokapot_paths = [
        os.path.join(mokapot_model_folder, f"mokapot_model_{i}.pkl") for i in range(3)
    ]

    calibration_folder = "{}/{}/calibration".format(
        output_folder,
        run_name
    )

    os.makedirs(os.path.dirname(save_feature_path), exist_ok=True)
    os.makedirs(os.path.dirname(save_psm_path), exist_ok=True)
    os.makedirs(calibration_folder, exist_ok=True)
    os.makedirs(mokapot_model_folder, exist_ok=True)

    return {
        'feature_path': save_feature_path,
        'psm_path': save_psm_path,
        'mokapot_model_folder': mokapot_model_folder,
        'mokapot_model_paths': save_mokapot_paths,
        'calibration_folder': calibration_folder
    }

def filter_spectra_by_ids(psm_list, spectrum_id_list):
    n_before = len(psm_list)
    logging.info('Filtering spectra based on ground-truth labels...')

    spectrum_ids = psm_list['spectrum_id']

    psm_list = psm_list[
        np.in1d(
            spectrum_ids,
            spectrum_id_list
        )
    ]

    n_after = len(psm_list)
    logging.info('Filtered {} PSMs from spectrum_id_list. {} PSMs remaining.'.format(
        n_before-n_after,
        n_after
    ))
    return psm_list