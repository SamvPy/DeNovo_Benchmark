import argparse
from glob import glob

from denovo_utils.parsers import DenovoEngineConverter
from denovo_utils.parsers.constants import EXTENSIONS
from denovo_utils.rescoring import DeNovoRescorer
from denovo_utils.rescoring.rescore import (
    load_configuration,
    load_rescorer,
    save_object,
    already_trained
)
import pickle

from ms2rescore import parse_configurations
import copy
import os
from glob import glob
import json

import logging

logger = logging.getLogger(__name__)


def parse_paths(
        ground_truth_folder,
        mgf_folder,
        denovo_folder,
        engines,
        *args,
        **kwargs
):
    
    mgf_psm_dict = {}
    for mgf_path in glob(os.path.join(mgf_folder, "*.mgf")):
        ground_truth_filetype = "sage"
        filename = os.path.basename(mgf_path).split(".")[0]
        
        psm_dict = {
            "psm_file": os.path.join(ground_truth_folder, filename + EXTENSIONS[ground_truth_filetype]),
            "spectrum_path": mgf_path
        }

        denovo_dict = {}
        for engine in engines:
            if engine in ["instanovoplus", "spectralis"]:
                for result_path in glob(
                    os.path.join(denovo_folder, engine)+f"/*/{filename}*"
                ):
                    base_denovo_engine = os.path.basename(os.path.dirname(result_path))
                    engine_label = f"{engine}__{base_denovo_engine}"
                    
                    denovo_dict[engine_label] = result_path
            else:
                denovo_dict[engine] = os.path.join(
                    denovo_folder, engine, filename + f".{engine}" + EXTENSIONS[engine]
                )
        psm_dict["denovo"] = denovo_dict
        mgf_psm_dict[filename] = psm_dict
    return mgf_psm_dict


def apply_pipeline(config, filename):

    logging.info(f"Process for {filename} started.")
    # Load and initialize the configurations
    logging.info("Reading configuration file")
    configuration = load_configuration(
        path_config_ms2rescore=config["path_ms2rescore_config"],
        psm_file=config["psm_file"],
        spectrum_path=config["spectrum_path"]
    )

    save_feature_path = "{}/{}/{}/{}.parquet".format(
        config["output_folder"],
        filename,
        "features",
        config["ground_truth_engine"]
    )
    save_psm_path = "{}/{}/{}/{}.pkl".format(
        config["output_folder"],
        filename,
        "psmlist",
        config["ground_truth_engine"]
    )

    mokapot_model_folder = "{}/{}/mokapot".format(
        config["output_folder"],
        filename
    )
    save_mokapot_paths = [
        os.path.join(mokapot_model_folder, f"mokapot_model_{i}.pkl") for i in range(3)
    ]

    # 0. Initialize rescoring object
    rescorer = DeNovoRescorer(configuration=configuration)

    # Check if models were already trained or not.
    # If so, it is pointless to rerun this part again.
    trained = already_trained(
        save_psm_path,
        save_feature_path,
        save_mokapot_paths
    )
    if trained:
        logging.info(f"Loading the mokapot models for {filename} from {mokapot_model_folder}")
        psm_list = load_rescorer(
            rescorer=rescorer,
            psm_path=save_psm_path,
            feature_path=save_feature_path,
            mokapot_folder=mokapot_model_folder
        )
    
    else:
        # 1. Load the ground-truth dataset
        logging.info(f"No mokapot models, rescored psm_list or fgen file found for {filename}")
        logging.info("Reading ground-truth data")
        parser = DenovoEngineConverter.select(config["ground_truth_engine"])
        psm_list = parser.parse(
            result_path=config["psm_file"],
            mgf_path=config["spectrum_path"]
        )
        psm_list = psm_list.get_rank1_psms()

        # 2. Preprocess the psm_list for rescoring
        psm_list = rescorer.preprocess_psm_list(psm_list)

        # 3. Set up the deeplc model by transfer learning
        logging.info("Retraining DeepLC")
        rescorer.retrain_deeplc(psm_list)

        # 4. Generate features for x psm_lists
        logging.info("Adding features for ground truth")
        rescorer.add_features(psm_list)

        # 5. Use the ground-truth for mokapot model training
        os.makedirs(mokapot_model_folder, exist_ok=True)
        logging.info("Training mokapot models")
        rescorer.train_mokapot_models(
            psm_list,
            save_folder=mokapot_model_folder
        )

        # 6. Apply the mokapot models for psm-scoring on the ground truth
        logging.info(f"Rescoring ground-truth. {len(psm_list)} psms")
        rescorer.rescore(psm_list)

        psm_features = psm_list.to_dataframe()[["spectrum_id", "run", "source", "rescoring_features"]]
        psm_features.to_parquet(save_feature_path)
        psm_list["rescoring_features"] = [{}]*len(psm_list)
        
        save_object(
            psm_list,
            save_psm_path,
            'parquet'
        )

    # 7. Apply the trained fgens and mokapot model on the denovo PSMS
    for engine, path in config["denovo"].items():
        # 7.1 Load the data
        # Add refiner token to the filename if required
        base_denovo = "base"
        if engine.startswith('spectralis') or engine.startswith("instanovoplus"):
            engine, base_denovo =  engine.split("__")
            save_feature_path = "{}/{}/{}/{}.{}.parquet".format(
                config["output_folder"],
                filename,
                "features",
                base_denovo,
                engine
            )
            save_psm_path = "{}/{}/{}/{}.{}.pkl".format(
                config["output_folder"],
                filename,
                "psmlist",
                base_denovo,
                engine
            )
        else:
            save_feature_path = "{}/{}/{}/{}.parquet".format(
                config["output_folder"],
                filename,
                "features",
                engine
            )
            save_psm_path = "{}/{}/{}/{}.pkl".format(
                config["output_folder"],
                filename,
                "psmlist",
                engine
            )

        if os.path.exists(save_psm_path):
            logging.info(f"Skipping {filename}:{engine}-{base_denovo} rescoring. {save_psm_path} exists.")
            continue

        
        logging.info(f"{engine}-{base_denovo} for {filename}")
        parser = DenovoEngineConverter.select(engine)
        psm_list = parser.parse(
            result_path=path,
            mgf_path=config["spectrum_path"]
        )
        
        # 7.2 Preprocess denovo dataset
        logging.info(f"Adding features for {engine}:{filename}")
        psm_list = rescorer.preprocess_psm_list(psm_list, denovo=True)
        rescorer.add_features(psm_list)

        # 7.3 Rescore de novo sequences
        logging.info(f"Rescoring for {engine}:{filename}. {len(psm_list)} psms")
        rescorer.rescore(psm_list, denovo=True)

        # 7.4 Split and save the PSMList in two distinct parts
        psm_features = psm_list.to_dataframe()[["spectrum_id", "run", "source", "rescoring_features"]]
        psm_features.to_parquet(save_feature_path)
        psm_list["rescoring_features"] = [{}]*len(psm_list)

        save_object(
            psm_list,
            save_psm_path,
            'pickle'
        )

def complete_config(config, args):
    config = copy.deepcopy(config)
    config["path_ms2rescore_config"] = args["path_ms2rescore_config"]
    config["ground_truth_engine"] = args["ground_truth_engine"]
    config["output_folder"] = args["output_folder"]
    return config

def main(args):
    
    mgf_psm_dict = parse_paths(**args)

    if args["filename"] in mgf_psm_dict.keys():
        # Only run pipeline for that file
        config = complete_config(
            config=mgf_psm_dict[args["filename"]],
            args=args
        )
        apply_pipeline(config, filename)

    for filename, paths in mgf_psm_dict.items():
        config = complete_config(
            config=paths,
            args=args
        )
        apply_pipeline(config, filename)



def load_config(config_path):
    with open(config_path, 'r') as f:
        return json.load(f)

if __name__ == "__main__":

    # change to config file
    parser = argparse.ArgumentParser(description="Run MS2Rescore de-novo pipeline.")
    parser.add_argument(
        "-c", "--config", default=None, help="Path to JSON configuration file"
    )
    parser.add_argument(
        "-m", "--mgf_folder", required=False, help="Folder to mgf files"
    )
    parser.add_argument(
        "-g", "--ground_truth_folder", required=False, help="Folder to ground_truth files"
    )
    parser.add_argument(
        "-d", "--denovo_folder", required=False
    )
    parser.add_argument(
        "-r", "--path_ms2rescore_config", required=False, help="Path to the ms2rescore pipeline config"
    )
    parser.add_argument(
        "-e", "--ground_truth_engine", required=False, help="The name of the ground truth engine"
    )
    parser.add_argument(
        "-a", "--engines", required=False, help="Denovo engine output to rescore"
    )
    parser.add_argument(
        "-o", "--output_folder", required=False, help="Output folder to store psm lists and mokapot models"
    )
    parser.add_argument(
        "-t", "--file_type", default="dir", help="dir or file?"
    )
    parser.add_argument(
        "-f", "--filename", default=None, help="Specify filename of only one file should be processed"
    )

    args = parser.parse_args()

    # Load config if provided
    if args.config:
        
        with open(args.config, 'r') as f:
            config_args = json.load(f)
        # Override command-line arguments with values from the config file if they aren't provided via the CLI
        for key, value in config_args.items():
            if getattr(args, key) is None:  # If the argument is not provided in the CLI, override with the config
                setattr(args, key, value)

    logging.info("Provided arguments: {}".format(vars(args)))
    main(vars(args))