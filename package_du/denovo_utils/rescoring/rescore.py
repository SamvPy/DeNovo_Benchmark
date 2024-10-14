# Largely based on ms2rescore codebase (core.py)
# Should allow to interface with the feature generators once more.
# SHould make it possible to store and load feature generators to and from a file.
# I dont know if this is even possible considering models will be need to stored aswell... 
# Especially for DeepLC...

import json
import logging
from multiprocessing import cpu_count
from typing import Dict, Optional
import copy

import os
import numpy as np
import psm_utils.io
from mokapot.dataset import LinearPsmDataset
from psm_utils import PSMList
import pandas as pd

from ms2rescore import exceptions
from ms2rescore.feature_generators import FEATURE_GENERATORS
from ms2rescore.parse_psms import parse_psms
from ms2rescore.parse_spectra import get_missing_values
from ms2rescore.report import generate
from ms2rescore.rescoring_engines import mokapot, percolator
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

class FGens:
    def __init__(self, configuration):
        self.config = configuration

    def initialize_feature_generators(self):
        
        self.feature_generators = {}
        for fgen_name, fgen_config in self.config["feature_generators"].items():
            if fgen_name == "maxquant":
                continue
            conf = self.config.copy()
            conf.update(fgen_config)
            fgen = FEATURE_GENERATORS[fgen_name](**conf)

            self.feature_generators[fgen_name] = fgen
    

class DeNovoRescorer:
    def __init__(self, configuration):
        self.configuration = configuration
        self.fgens = FGens(configuration=configuration["ms2rescore"])
        self.fgens.initialize_feature_generators()


    def preprocess_psm_list(
            self, psm_list, denovo=False
    ):
        config = self.configuration["ms2rescore"]
        psm_list = parse_psms(config, psm_list, recalculate_qvalues=not denovo)
        psm_list = _fill_missing_precursor_info(psm_list, config)
        if denovo:
            psm_list["retention_time"] = psm_list["retention_time"]/60

            # MS2PIP fill fail if the sequence is shorter than 4 amino acids
            mask = np.array([
                len(psm.peptidoform)>=4 for psm in psm_list
            ])
            dropped = np.sum(~mask)
            if dropped > 0:
                logging.warning(f"Dropped {dropped} spectra due to peptides with length <4")

            psm_list = psm_list[mask]
        return psm_list

    def retrain_deeplc(self, psm_list):
        # The model is saved in a temp folder so should edit here to save the deeplc model directly.
        self.fgens.feature_generators["deeplc"].retrain_deeplc(psm_list)

    def add_features(self, psm_list):
        if self.fgens.feature_generators["deeplc"].deeplc_predictor is None:
            raise Exception("First retrain deeplc!")
            
        config = self.configuration["ms2rescore"]
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

    def filter_psm_list(self, psm_list, denovo=False):
        psm_list_rescoring = copy.deepcopy(psm_list)
        psm_list_rescoring["rescoring_features"] = np.array(list(pd.DataFrame(psm_list_rescoring["rescoring_features"].tolist())[
            self.configuration["ms2rescore"]["rescoring_features"]
        ].to_dict(orient="index").values()))

        # Otherwise cannot convert to a mokapot PSMDataset
        if denovo:
            half_targets = [True]*(int(len(psm_list_rescoring)/2))
            half_decoys = [False]*(len(psm_list_rescoring)-len(half_targets))
            psm_list_rescoring["is_decoy"] = half_targets+half_decoys

        return psm_list_rescoring

    def train_mokapot_models(self, psm_list, save_folder=None):
        """
        Use psm_list and configuration to train mokapot models. They will be stored as attributes
        and are optionally saved to a given folder (save_folder)
        """
        # Bunch of config code
        _set_log_levels()

        # Filter to features specified
        psm_list_rescoring = self.filter_psm_list(psm_list)

        config = self.configuration["ms2rescore"]
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
                lin_psm_data, model=PercolatorModel(train_fdr=config["rescoring_engine"]["mokapot"]["train_fdr"]), rng=8
            )
        except RuntimeError as e:
            raise exceptions.RescoringError("Mokapot could not be run. Please check the input data.") from e

        if save_folder is not None:
            for i, model in enumerate(models):
                model.save(
                    os.path.join(save_folder, f"mokapot_model_{i}.pkl")
                )
        self.mokapot_models = models


    def rescore(self, psm_list, denovo=False):
        """
        Rescore using trained mokapot models.
        """
        config = self.configuration["ms2rescore"]
        # Rescore PSMs
        _set_log_levels()

        psm_list_rescoring = self.filter_psm_list(psm_list, denovo=denovo)
        lin_psm_data = convert_psm_list(psm_list_rescoring)

        # The brew function performs a calibration between the fold based on the scores outputted.
        # I really dont want this because if i use the same models but a different psm_list, this
        # is scored with a differing scoring function!!!

        scores = _predict(
            lin_psm_data, models=self.mokapot_models
        )

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

        else:
            for psm in psm_list:
                psm["provenance_data"].update(
                    {f"score_{psm.source}": psm["score"]}
                )
            psm_list["score"] = scores

    
def _predict(dset, models):
    """
    Return the new scores for the dataset

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