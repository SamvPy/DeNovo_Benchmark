import argparse
from glob import glob

from denovo_utils.parsers import DenovoEngineConverter
from denovo_utils.parsers.constants import EXTENSIONS
from denovo_utils.rescoring import DeNovoRescorer
from denovo_utils.io.read import (
    load_configuration,
    load_psmlist_and_features
)
from denovo_utils.io.save import save_pickle, save_features, save_psmlist

from denovo_utils.rescoring.utils import (
    already_trained,
    parse_config_paths,
    setup_paths,
    filter_spectra_by_ids
)
import gc

from ms2rescore import parse_configurations
import copy
import os
from glob import glob
import json

import logging

logger = logging.getLogger(__name__)


def apply_pipeline(config, filename, args=None):

    logging.info(f"Process for {filename} started.")
    # Load and initialize the configurations
    logging.info("Reading configuration file")
    configuration = load_configuration(
        path_config_ms2rescore=config["path_ms2rescore_config"],
        psm_file=config["psm_file"],
        spectrum_path=config["spectrum_path"]
    )

    config['save_paths'] = setup_paths(
        output_folder=config['output_folder'],
        run_name=os.path.basename(config['psm_file']).split('.')[0],
        filename='ground_truth'
    )

    # 0. Initialize rescoring object
    rescorer = DeNovoRescorer(configuration=configuration)

    # Check if models were already trained or not.
    # If so, it is pointless to rerun this part again.
    trained = already_trained(
        config['save_paths']['psm_path'],
        config['save_paths']['feature_path'],
        config['save_paths']['mokapot_model_paths']
    )
    if trained:
        logging.info(f"Found result paths for {filename} ground truth")
        rescorer.load(
            calibration_folder=config['save_paths']['calibration_folder'],
            mokapot_folder=config['save_paths']['mokapot_model_folder']
        )

        psm_list = load_psmlist_and_features(
            psm_path=config['save_paths']['psm_path'],
            feature_path=config['save_paths']['feature_path']
        )
        psm_list = psm_list.get_rank1_psms()

        if len(args["regenerate"]) > 0:
            logging.info(f"Regenerating {args['regenerate']} features.")
            rescorer.add_features(psm_list=psm_list, feature_types=args['regenerate'])

        if args['retrain_mokapot']:
            logging.info("Training mokapot models")
            logging.info(f"Retrain mokapot and rescore ground-truth. {len(psm_list)} psms")
            rescorer.train_mokapot_models(psm_list, save_folder=config['save_paths']['mokapot_model_folder'])
            rescorer.rescore(psm_list)

            save_features(
                psm_list,
                config['save_paths']['feature_path']
            )            
            save_psmlist(
                psm_list,
                config['save_paths']['psm_path']
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

        # 3. Calibrate feature generators
        logging.info("Calibrating feature generators...")
        rescorer.calibrate_fgens(
            psm_list=psm_list,
            calibration_folder=config['save_paths']['calibration_folder'],
            from_cache=False
        )

        # 4. Generate features for x psm_lists
        logging.info("Adding features for ground truth")
        rescorer.add_features(psm_list)

        # 5. Use the ground-truth for mokapot model training
        logging.info("Training mokapot models")
        rescorer.train_mokapot_models(
            psm_list,
            save_folder=config['save_paths']['mokapot_model_folder']
        )

        # 6. Apply the mokapot models for psm-scoring on the ground truth
        logging.info(f"Rescoring ground-truth. {len(psm_list)} psms")
        rescorer.rescore(psm_list)

        save_features(
            psm_list,
            config['save_paths']['feature_path']
        )            
        save_psmlist(
            psm_list,
            config['save_paths']['psm_path']
        )

    # Keep track of the spectrum_ids for the ground-truth dataset.
    spectrum_ids_gt = psm_list['spectrum_id']

    # 7. Apply the trained fgens and mokapot model on the denovo PSMS
    for engine, path in config["denovo"].items():
        # 7.1 Load the data
        # Add refiner token to the filename if required
        base_denovo = "base"
        if engine.startswith('spectralis') or engine.startswith("instanovoplus"):
            engine, base_denovo = engine.split("__")
            save_feature_path = "{}/{}/{}/{}.{}.parquet".format(
                config["output_folder"],
                filename,
                "features",
                base_denovo,
                engine
            )
            save_psm_path = "{}/{}/{}/{}.{}.parquet".format(
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
            save_psm_path = "{}/{}/{}/{}.parquet".format(
                config["output_folder"],
                filename,
                "psmlist",
                engine
            )

        # Load the data depending on how its specified in the config
        # This is more for dev purposes when feature generators are added/changed or when rescoring should be redone
        logging.info(f"{engine}-{base_denovo} for {filename}")
        if os.path.exists(save_psm_path) and os.path.exists(save_feature_path):
            
            if len(args["regenerate"]) > 0:
                logging.info(f"Regenerating {args['regenerate']} features.")
                psm_list = load_psmlist_and_features(psm_path=save_psm_path, feature_path=save_feature_path)
                rescorer.add_features(psm_list=psm_list, feature_types=args['regenerate'])
            elif not args['retrain_mokapot']:
                logging.info(f"Skipping {filename}:{engine}-{base_denovo} rescoring. {save_psm_path} exists.")
                continue
            else:
                psm_list = load_psmlist_and_features(psm_path=save_psm_path, feature_path=save_feature_path)
            logging.info(f"Loaded {filename} from {os.path.dirname(save_psm_path)}")

        # Run this part if no rescored psmlist and features for the denovo psms are stored in save paths.
        else:
            logging.info(f"Loading {filename} from raw search results.")
            parser = DenovoEngineConverter.select(engine)
            psm_list = parser.parse(
                result_path=path,
                mgf_path=config["spectrum_path"]
            )
            
            # Filter out spectra which cannot be compared to ground-truth. Saves computation time
            if args['only_gt_spectra']:
                psm_list = filter_spectra_by_ids(
                    psm_list=psm_list,
                    spectrum_id_list=spectrum_ids_gt
                )

            # 7.2 Preprocess denovo dataset
            if not args['batches']:
                logging.info(f"Adding features for {engine}:{filename}")
                psm_list = rescorer.preprocess_psm_list(psm_list, denovo=True)
                rescorer.add_features(psm_list)

        if args['batches']:
            rescorer.batched_fgen_and_rescore(
                psm_list=psm_list,
                n_psms=10000,
                denovo=True,
                save_folder=os.path.dirname(os.path.dirname(save_feature_path)),
                filename=os.path.basename(save_feature_path.split('.')[0])
            )

        else:
            # 7.3 Rescore de novo sequences
            logging.info(f"Rescoring for {engine}:{filename}. {len(psm_list)} psms")
            rescorer.rescore(psm_list, denovo=True)

            # 7.4 Split and save the PSMList in two distinct parts
            save_features(
                psm_list,
                save_feature_path
            )            
            save_psmlist(
                psm_list,
                save_psm_path
            )
        gc.collect()

def complete_config(config, args):
    config = copy.deepcopy(config)
    config["path_ms2rescore_config"] = args["path_ms2rescore_config"]
    config["ground_truth_engine"] = args["ground_truth_engine"]
    config["output_folder"] = args["output_folder"]
    return config

def main(args):
    
    mgf_psm_dict = parse_config_paths(**args)

    if args["filename"] in mgf_psm_dict.keys():
        # Only run pipeline for that file
        config = complete_config(
            config=mgf_psm_dict[args["filename"]],
            args=args
        )
        apply_pipeline(config, filename, args=args)

    for filename, paths in mgf_psm_dict.items():
        if filename not in args["filenames"]:
            logging.info(f"Skipping processing for {filename}.")
            continue

        config = complete_config(
            config=paths,
            args=args
        )
        apply_pipeline(config, filename, args=args)


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
    parser.add_argument(
        "-b", "--batches", default=False, help='Whether to perform rescoring of de novo results in batches.'
    )

    args = parser.parse_args()

    # Load config if provided
    args = vars(args)

    if "config" in args.keys():
        
        with open(args['config'], 'r') as f:
            config_args = json.load(f)
        # Override command-line arguments with values from the config file if they aren't provided via the CLI
        for key, value in config_args.items():
            # try:
            #     if args[key] is None:  # If the argument is not provided in the CLI, override with the config
            #         args[key] = value
            # except:
            #     args[key] = value
            args[key] = value

    logging.info("Provided arguments: {}".format(args))
    main(args)