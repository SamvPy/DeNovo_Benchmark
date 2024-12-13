import logging
from typing import Optional

import shutil
import os
from psm_utils import PSMList

from ..io.read import load_pickle, load_deeplc
from ..io.save import save_pickle

from tensorflow.keras.models import load_model

from ms2rescore.feature_generators import FEATURE_GENERATORS

logger = logging.getLogger(__name__)

class FGens:
    """
    Store and initialize feature generators in a central object.

    Attributes
    ----------
    config: dict
        Contains settings for each feature generator. (see ms2rescore config for reference)
    feature_generators: dict
        keys: feature generator name
        values: feature generator object
    
    Method
    ------
    initialize_feature_generators:
        Initialize feature generators with provided configuration file.
    """
    def __init__(self, configuration):
        """
        Initialize feature generator group object.
        
        Parameter
        ----------
        config: dict
            Contains settings for each feature generator. (see ms2rescore config for reference)
        """
        self.config = configuration
        self.deeplc_calibration_name = "deeplc_calibration.pkl"
        self.im2deep_calibration_name = "im2deep_calibration.pkl"

    def initialize_feature_generators(self):
        """
        Initialize feature generators according to provided configuration file.

        Parameter
        ---------
        feature_generators: dict
            keys: feature generator name
            values: feature generator object
        """
        self.feature_generators = {}
        for fgen_name, fgen_config in self.config["feature_generators"].items():
            if fgen_name == "maxquant":
                continue
            conf = self.config.copy()
            conf.update(fgen_config)
            fgen = FEATURE_GENERATORS[fgen_name](**conf)

            self.feature_generators[fgen_name] = fgen

    def calibrate(
            self,
            calibration_folder: str,
            from_cache: bool,
            psm_list: Optional[PSMList]=None,
    ) -> None:
        """
        Calibrate the feature generators using psm_list.

        When calibrating, everything required for calibration (models and calibration files)
        are stored to the calibration folder. Alternatively, they are loaded from this location
        when using caching.

        Currently only im2deep and deeplc require calibration.

        Parameters
        ----------
        psm_list: PSMList
            The psm_list used for calibration
        calibration_folder: str
            The calibration related files/models will be stored in this folder
        from_cache: bool
            If True, will calibrate fgens by loading the relevant files. Else will calibrate from scratch
            and store the calibration data in calibration_folder.
        """
        if psm_list is None and not from_cache:
            raise Exception("Cannot calibrate fgens from scratch while not having a psm_list.")

        if 'deeplc' in self.feature_generators.keys():
            # Specify path for storage
            deeplc_calibration_path=os.path.join(
                calibration_folder,
                self.deeplc_calibration_name
            )

            # Load or create/save
            if from_cache:
                self.feature_generators['deeplc'] = load_deeplc(
                    model_path_folder=calibration_folder,
                    calibration_path=deeplc_calibration_path
                )
            else:
                self._retrain_deeplc(
                    psm_list,
                    save_model_folder=calibration_folder
                )
                # Reload the DeepLC model to prevent it from 
                # using the model from the temporary directory.
                deeplc_fgen = load_deeplc(
                    model_path_folder=calibration_folder,
                    calibration_path=deeplc_calibration_path
                )
                self.fgens.feature_generators['deeplc'] = deeplc_fgen
        
        if 'im2deep' in self.feature_generators.keys():
            # Specify path for storage
            im2deep_calibration_path = os.path.join(
                calibration_folder,
                self.im2deep_calibration_name
            )
            
            # Load or create/save
            if from_cache:
                im2deep_calibration_df = load_pickle(im2deep_calibration_path)
            else:
                im2deep_calibration_df = self.feature_generators["im2deep"].make_calibration_df(
                    psm_list=psm_list
                )
                save_pickle(im2deep_calibration_df, im2deep_calibration_path)
            self.feature_generators['im2deep'].cal_psm_df = im2deep_calibration_df

    def _retrain_deeplc(
            self,
            psm_list: Optional[PSMList]=None,
            save_model_folder: str=""
        ) -> None:
        """
        Fine tune a deeplc model using transfer learning and calibration.
        """
        # Perform retraining (transfer learning) of the deeplc model using psm_list
        # The model is saved in a temp folder. Will additionally save to save_model_folder.
        self.feature_generators["deeplc"].retrain_deeplc(psm_list)
        
        # Save the deeplc model checkpoint
        deeplc_fgen = self.fgens.feature_generators["deeplc"]
        for m in deeplc_fgen.selected_model:
            new_model_path = os.path.join(save_model_folder, os.path.basename(m))
            
            loaded_model = load_model(m)
            loaded_model.save(new_model_path)

            # Delete the temporary directory (Hopes to prevent broken pipes)
            shutil.rmtree(os.path.dirname(m))  # Manually remove a temporary folder
            logging.info(f"Saved DeepLC model to {new_model_path}")

        # Create calibration parameters and save it
        calibration_params = {
            "calibrate_dict": {new_model_path: deeplc_fgen.deeplc_predictor.calibrate_dict[m]},
            "calibrate_min": {new_model_path: deeplc_fgen.deeplc_predictor.calibrate_min[m]},
            "calibrate_max": {new_model_path: deeplc_fgen.deeplc_predictor.calibrate_max[m]}
        }
        calibration_path = os.path.join(save_model_folder, self.deeplc_calibration_name)
        save_pickle(calibration_params, calibration_path)
        logging.info(f"Saved DeepLC calibration to {calibration_path}")