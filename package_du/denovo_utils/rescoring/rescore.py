# Largely based on ms2rescore codebase (core.py)
# Should allow to interface with the feature generators once more.
# SHould make it possible to store and load feature generators to and from a file.
# I dont know if this is even possible considering models will be need to stored aswell... 
# Especially for DeepLC...

import json
import logging
from typing import Dict, Optional
import copy
from glob import glob
import os

import numpy as np
import psm_utils.io
from psm_utils import PSMList
import pandas as pd
from ..io.save import save_mokapot_models, save_pickle, save_features, save_psmlist
from ..io.read import load_mokapot
from .fgens import FGens

from ms2rescore import exceptions
from ms2rescore.parse_psms import parse_psms
from ms2rescore.rescoring_engines.mokapot import add_peptide_confidence, add_psm_confidence

from ms2rescore.core import (
    _fill_missing_precursor_info,
    _write_feature_names,
    _fix_constant_pep,
    _filter_by_rank,
    _calculate_confidence

)
from ms2rescore.rescoring_engines.mokapot import (
    _set_log_levels,
    convert_psm_list,
    add_peptide_confidence,
    add_psm_confidence
)

from mokapot import brew
from mokapot.model import PercolatorModel

logger = logging.getLogger(__name__)


class DeNovoRescorer:
    """
    Rescoring object.

    Encapsulates all data and functionalities required to rescore a de novo psmlist.
    """
    def __init__(self, configuration):
        """
        Initialize DeNovoRescorer by initializing the feature generators.

        Parameter
        ---------
        configuration: dict
            Contains settings for each feature generator. (see ms2rescore config for reference)
        """
        self.config = configuration
        self.fgens = FGens(configuration=configuration["ms2rescore"])
        self.fgens.initialize_feature_generators()

    def load(
            self,
            calibration_folder,
            mokapot_folder=None
    ):
        """
        Load all essential data for rescoring from a rescorer created in the past.

        Required data for a reproducible rescorer are:
        - calibrated feature generators (to ensure reproducible feature generation)
        - mokapot models (to ensure same rescoring function)

        Parameters
        ----------
        calibration_folder: str
            The folder where feature generation-related files are stored.
        mokapot_folder: str
            The folder where mokapot models are stored.
        """
        self.fgens.calibrate(
            calibration_folder=calibration_folder,
            from_cache=True
        )

        if mokapot_folder is not None:
            self.mokapot_models = load_mokapot(mokapot_folder)

    def calibrate_fgens(
            self,
            psm_list,
            calibration_folder,
            from_cache: bool=False
        ) -> None:
        """
        Calibrate the feature generators.
        """
        # Preprocess the psm_list accordingly. Denovo results should not be used for calibration
        # THOUGHT: Maybe allow calibrating on denovo results ? This would allow ambiguity estimation
        # without requiring a ground-truth dataset ? Similar to NOKOI?
        psm_list = self.preprocess_psm_list(psm_list=psm_list, denovo=False)

        # Calibrate the feature generators
        self.fgens.calibrate(
            psm_list=psm_list,
            calibration_folder=calibration_folder,
            from_cache=from_cache
        )

    def preprocess_psm_list(
            self, psm_list: PSMList, denovo: bool=False
    ) -> PSMList:
        """
        Preprocess psm_list prior to rescoring.

        This process mainly entails the following:
        - Filter psms with ranks lower than configured in max_psm_rank_input
        - Remove peptidoforms with invalid amino acids
        - Recalculate q-values if not denovo psmlist and if q-values are NOT present
        - Move score (qvalue, pep and rank) to provenance data
        - Parse modifications as configured in modification_mapping and fixed_modifications
        - For denovo: Filter out peptidoforms with length < 4 and add random is_decoy labels

        Parameters
        ----------
        psm_list: PSMList
            PSMList to preprocess.
        denovo: bool
            Whether the PSMList contains de novo psms and is thus not used as ground-truth.
            If True, will not recalculate q-values and do a filtering step on peptidoform length
            and add is_decoy labels randomly, as this is required when using mokapot.

        Return
        ------
        PSMList
            Preprocessed PSMList.
        """
        config = self.config["ms2rescore"]

        # max_psm_rank_input is important in the config when 
        # rescoring multiple de novo hits of multiple ranks.
        psm_list = parse_psms(config, psm_list, recalculate_qvalues=not denovo)
        psm_list = _fill_missing_precursor_info(psm_list, config)
        
        # De novo specific preprocessing
        if denovo:
            # Filter psms based on length
            # MS2PIP fill fail if the sequence is shorter than 4 amino acids
            mask = np.array([
                len(psm.peptidoform)>=4 for psm in psm_list
            ])
            dropped = np.sum(~mask)
            if dropped > 0:
                logging.warning(f"Dropped {dropped} spectra due to peptides with length <4")
            psm_list = psm_list[mask]

            # Add is_decoy labels randomly
            # Otherwise parsing to a PSMDataset in mokapot is not possible.
            half_targets = [True]*(int(len(psm_list)/2))
            half_decoys = [False]*(len(psm_list)-len(half_targets))
            psm_list["is_decoy"] = half_targets+half_decoys

        return psm_list


    def add_features(self, psm_list: PSMList, feature_types: list=[]) -> None:
        """
        Add features to the psm_list using the feature generators loaded in the rescorer object.

        Parameters
        ----------
        psm_list: PSMList
            Add features to psm_list in the 'rescoring_feature' field.
        feature_types: list
            Specifies which features to add. Corresponds to the name of the feature generator.
            Defaults to all feature generators loaded in the rescorer.
            To select only a given set of features for rescoring, define this in ms2rescore config
            as this is handled elsewhere.
        """
        # Default feature types to all feature generators in the fgen object
        if feature_types == []:
            feature_types = list(self.fgens.feature_generators.keys())

        if self.fgens.feature_generators["deeplc"].deeplc_predictor is None:
            raise Exception("First retrain deeplc!")
            
        config = self.config["ms2rescore"]
        output_file_root = config["output_path"]

        # Define feature names; get existing feature names from PSM file
        feature_names = dict()
        psm_list_feature_names = {
            feature_name
            for psm_list_features in psm_list["rescoring_features"]
            for feature_name in psm_list_features.keys()
        }
        feature_names["psm_file"] = psm_list_feature_names
        logger.debug(
            f"PSMs already contain the following rescoring features: {psm_list_feature_names}"
        )

        # Add rescoring features
        for fgen_name, fgen_config in config["feature_generators"].items():
        # TODO: Handle this somewhere else, more generally?
            if fgen_name == "maxquant" and not (psm_list["source"] == "msms").all():
                logger.warning(
                    "MaxQuant feature generator requires PSMs from a MaxQuant msms.txt file. Skipping "
                    "this feature generator."
                )
                continue
            if fgen_name not in feature_types:
                logging.info(f"Skipping feature generators of type {fgen_name}")
                continue

            logger.info(f"Generating {fgen_name} features")

            fgen = self.fgens.feature_generators[fgen_name]
            fgen.add_features(psm_list)

            logger.debug(f"Adding features from {fgen_name}: {set(fgen.feature_names)}")
            feature_names[fgen_name] = set(fgen.feature_names)

        # Filter out psms that do not have all added features
        all_feature_names = {f for fgen in feature_names.values() for f in fgen}
        psms_with_features = [
            (set(psm.rescoring_features.keys()) == all_feature_names) for psm in psm_list
        ]

        if psms_with_features.count(False) > 0:
            removed_psms = psm_list[[not psm for psm in psms_with_features]]
            missing_features = {
                feature_name
                for psm in removed_psms
                for feature_name in all_feature_names - set(psm.rescoring_features.keys())
            }
            logger.warning(
                f"Removed {psms_with_features.count(False)} PSMs that were missing one or more "
                f"rescoring feature(s), {missing_features}."
            )
        psm_list = psm_list[psms_with_features]  

        # Write feature names to file
        _write_feature_names(feature_names, output_file_root)

        if config["rename_to_usi"]:
            logging.debug(f"Creating USIs for {len(psm_list)} PSMs")
            psm_list["spectrum_id"] = [psm.get_usi(as_url=False) for psm in psm_list]

        # If no rescoring engine is specified or DEBUG, write PSMs and features to PIN file
        if not config["rescoring_engine"] or config["log_level"] == "debug":
            logger.info(f"Writing added features to PIN file: {output_file_root}.psms.pin")
            psm_utils.io.write_file(
                psm_list,
                output_file_root + ".pin",
                filetype="percolator",
                feature_names=all_feature_names,
            )

        if not config["rescoring_engine"]:
            logger.info("No rescoring engine specified. Skipping rescoring.")
            return None

    def filter_psm_list(self, psm_list: PSMList) -> PSMList:
        """
        Filter the rescoring features in the psm_list to only those specified in the config.

        Config field used as filter: ['ms2rescore']['rescoring_features']

        Parameter
        ---------
        psm_list: PSMList
            psm_list for which features should be filtered. The filter is specified in loaded config.
        
        Return
        ------
        PSMList
            PSMList with filtered rescoring_features field.
        """
        psm_list_rescoring = copy.deepcopy(psm_list)
        psm_list_rescoring["rescoring_features"] = np.array(
            list(
                pd.DataFrame(
                    psm_list_rescoring["rescoring_features"].tolist() # Parse feature field
                )[
                    self.config["ms2rescore"]["rescoring_features"] # Filter features
                ].to_dict(orient="index").values()
            )
        )
        return psm_list_rescoring

    def train_mokapot_models(self, psm_list: PSMList, save_folder=None):
        """
        Use psm_list and configuration to train mokapot models. They will be stored as attributes
        and are optionally saved to a given folder (save_folder)
        """
        # Bunch of config code
        _set_log_levels()

        # Filter to features specified
        psm_list_rescoring = self.filter_psm_list(psm_list)

        config = self.config["ms2rescore"]
        if "fasta_file" not in config["rescoring_engine"]["mokapot"]:
            config["rescoring_engine"]["mokapot"]["fasta_file"] = config["fasta_file"]
        if "protein_kwargs" in config["rescoring_engine"]["mokapot"]:
            kwargs = config["rescoring_engine"]["mokapot"].pop("protein_kwargs")
        else:
            kwargs = dict()

        # Parsing psm_list to mokapot compatible format
        lin_psm_data = convert_psm_list(psm_list_rescoring)

        # Rescore
        logger.debug(f"Mokapot brew options: `{kwargs}`")
        try:
            _, models = brew(
                lin_psm_data,
                model=PercolatorModel(
                    train_fdr=config["rescoring_engine"]["mokapot"]["train_fdr"]
                ),
                rng=8
            )
        except RuntimeError as e:
            raise exceptions.RescoringError("Mokapot could not be run. Please check the input data.") from e

        save_mokapot_models(models, save_folder)
        save_pickle(lin_psm_data.features.columns.tolist(), os.path.join(save_folder, 'feature_names.pkl'))
        self.mokapot_models = models


    def rescore(self, psm_list, denovo=False):
        """
        Rescore using trained mokapot models.
        """
        config = self.config["ms2rescore"]
        # Rescore PSMs
        _set_log_levels()

        psm_list_rescoring = self.filter_psm_list(psm_list)
        lin_psm_data = convert_psm_list(psm_list_rescoring)

        # The brew function performs a calibration between the fold based on the scores outputted.
        # I really dont want this because if i use the same models but a different psm_list, this
        # is scored with a differing scoring function!!!
        scores = _predict(
            lin_psm_data, models=self.mokapot_models
        )

        # When ground-truth, the previous confidence estimates will be updated.
        if not denovo:
            confidence = lin_psm_data.assign_confidence(
                scores=scores
            )
            add_psm_confidence(psm_list, confidence)
            add_peptide_confidence(psm_list, confidence)
            # Post-rescoring processing
            if all(psm_list["pep"] == 1.0):
                psm_list = _fix_constant_pep(psm_list)
            psm_list = _filter_by_rank(psm_list, config["max_psm_rank_output"], False)
            psm_list = _calculate_confidence(psm_list)

        # For denovo psm_lists, will add a single score to the provenance_data dictionary.
        else:
            for psm, ms2rescore_score in zip(psm_list, scores):
                psm["provenance_data"].update(
                    {f"score_ms2rescore": str(ms2rescore_score)}
                )

    def batched_fgen_and_rescore(self, psm_list, n_psms=10000, denovo=True, save_folder="", filename='filename'):
        # TODO: Implement a batched version of add_features and rescore
        # This should include a concatenation for saving psm_lists and features.
        os.makedirs(os.path.dirname(
            os.path.join(save_folder, 'features')
        ), exist_ok=True)
        os.makedirs(os.path.dirname(
            os.path.join(save_folder, "psmlist")
        ), exist_ok=True)

        counter = 0
        lower = 0
        upper = min(len(psm_list), n_psms)
        psm_list = self.preprocess_psm_list(psm_list, denovo=True)
        

        while True:
            # Select subset
            psm_subset = psm_list[lower: upper]

            # Preprocess again to make sure some decoy or targets
            # are present
            psm_subset = self.preprocess_psm_list(psm_subset, denovo=True)

            # Add features and rescore
            self.add_features(psm_subset)
            self.rescore(psm_subset, denovo=True)

            # Save
            save_features(
                psm_list=psm_subset,
                save_path=os.path.join(
                    save_folder,
                    "features",
                    f"{filename}_{counter}.parquet"
                )
            )
            save_psmlist(
                psm_list=psm_subset,
                save_path=os.path.join(
                    save_folder,
                    "psmlist",
                    f"{filename}_{counter}.parquet"
                )
            )

            # Break if the upper signified the end.
            if upper == len(psm_list):
                break
            lower += n_psms
            upper += n_psms
            upper = min(upper, len(psm_list))
            counter += 1


def _predict(dset, models):
    """
    Return the new scores for the dataset.

    Equals the average of the predictions from all models.

    Parameters
    ----------
    dset : PsmDataset
        The dataset to rescore
    models : list of Model
        The models for each dataset and whether it
        was reset or not.

    Returns
    -------
    numpy.ndarray
        A :py:class:`numpy.ndarray` containing the new scores.
    """
    test_set = copy.copy(dset)
    scores = []

    for model in models:
        scores.append(model.predict(test_set))
    
    return np.array(scores).mean(axis=0)
