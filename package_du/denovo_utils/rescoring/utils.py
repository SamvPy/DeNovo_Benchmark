import os
from glob import glob
import logging

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
            logger.warning(f'Multiple files matching {psm_file_path}. Taking {paths[0]}.')
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