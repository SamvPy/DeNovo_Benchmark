from psm_utils import PSMList, PSM
from psm_utils.io import read_file
import os
from ..constants import ENGINES, ENGINES_MAPPING
import pandas as pd
from pyteomics import mgf
from typing import Union
import numpy as np
from .base import DenovoEngineConverter
from tqdm import tqdm

tqdm.pandas()

class SpectralisParser():
    def __init__(self, mgf_path: str, results_dir: str):
        """
        Spectralis output parser towards PSMList.

        The class structure makes it possible to add results from disparate
        spectralis outputs to a single psmlist, making it easier to analyse.

        Parameters
        ----------
        mgf_path: str
            Path towards the mgf file analysed with spectralis.
        result_dir: str
            Path towards the folder storing all results from disparate de novo search engines

        Example
        -------
        ```python
        parser_spectralis = SpectralisParser(
                mgf_path=".../.mgf",
                results_dir=".../results"
            )
        parser_spectralis.add_xtandem(result_path="...")
        parser_spectralis.parse(path_spectralis="...")

        # Check which results were added
        parser_spectralis.added_results

        # Access each psmlist separately in a dictionary
        psmlist_dict = parser_spectralis.psmlists

        # Access a concatenated psmlist in dataframe format
        df = parser_spectralis.dataframe
        ```
        """
        self.mgf_path = mgf_path
        self.filename = os.path.basename(mgf_path).split(".")[0]
        self.results_dir = results_dir

        self.added_xtandem = False
        self.added_results = {
            engine: False for engine in ENGINES
        }
        self.psmlists = {}
    
    def add_xtandem(self, result_path):
        self._check_filename(result_path)
        self.psmlists["xtandem"] = read_file(
            result_path,
            filetype="percolator"
        )
        self.added_xtandem = True

    def parse(self, path_spectralis, overwrite=False) -> Union[PSMList, None]:
        columns = [
            "Spectralis_score",
            "scans",
            "filename"
        ]

        self._check_filename(path_spectralis)
        filename_spectralis = os.path.basename(path_spectralis).split(".")[0]
        
        spectralis_out = pd.read_csv(
            path_spectralis
        )
        
        spectralis_out_engines = spectralis_out.apply(
            reformat_engine_column, axis=1
        )
        scan_engine_df = pd.DataFrame(spectralis_out_engines.tolist(), columns=["scans", "engine"])
        spectralis_out_processed = pd.concat(
            [
                spectralis_out["Spectralis_score"],
                scan_engine_df["scans"],
                pd.DataFrame(scan_engine_df["engine"].tolist()),
            ], axis=1
        )

        spectralis_out_processed, loaded_engines = self._clean_dataframe(spectralis_out_processed)
        spectralis_out_processed["filename"] = filename_spectralis

        for engine in loaded_engines:
            if not overwrite and self.added_results[engine]:
                print(f"Already loaded results for {engine}! Skipping...")
                continue
                
            parser = DenovoEngineConverter.select(
                ENGINES_MAPPING[engine]
            )
            psm_df = parser.parse(
                result_path=os.path.join(
                    self.results_dir, ENGINES_MAPPING[engine], self.filename
                ),
                mgf_path=self.mgf_path
            ).to_dataframe()

            spectralis_engine_results = spectralis_out_processed.loc[
                spectralis_out_processed[engine],
                ["scans", "Spectralis_score"]
            ].set_index("scans")

            _ = psm_df.progress_apply(
                lambda x: self._add_spectralis_score(
                    x,
                    spectralis_engine_results
                ),
                axis=1
            )

            self.psmlists[ENGINES_MAPPING[engine]] = PSMList(
                psm_list = [PSM(**x) for x in psm_df.to_dict("records")]
            )
            self.added_results[engine] = True

    @property
    def psmlist(self) -> PSMList:
        psmlist_all = PSMList(psm_list=[])
        for psmlist in self.psmlists.values():
            psmlist_all += psmlist
        return psmlist_all
    
    @property
    def dataframe(self) -> pd.DataFrame:
        return self.psmlist.to_dataframe()
    
    def return_results(self, return_type: Union[PSMList, pd.DataFrame], engines: list):
        if isinstance(return_type, PSMList):
            pass

        elif isinstance(return_type, pd.DataFrame):
            pass

        else:
            raise TypeError("Provide the correct type in return_type! " \
                            f"Passed {type(return_type)}, only pd.DataFrame and PSMList are supported.")


    def _check_filename(self, p):
        p_filename = os.path.basename(p).split(".")[0]
        if self.filename not in p_filename:
            print(f"Warning! Are you sure you loaded the correct result file?")
            print(f"Loaded {p_filename} while mgf is {self.filename}.")

    def _clean_dataframe(self, df: pd.DataFrame):
        """
        Delete columns (de novo search engine column) which have no entries.
        """
        engine_empty = np.array(ENGINES)[
            (df.loc[:, ENGINES].sum()==0).tolist()
        ].tolist()

        loaded_engines = []
        for engine in ENGINES:
            if engine not in engine_empty:
                loaded_engines.append(engine)

        for engine in engine_empty:
            _ = df.pop(engine)
        return df, loaded_engines

    def _add_spectralis_score(self, row, spectralis_score_df):
        score = spectralis_score_df.loc[row["spectrum_id"], "Spectralis_score"]

        #TODO: If I want to include sequences not scored, this should return None
        try:
            row["rescoring_features"]["spectralis_score"] = float(score)
        except KeyError:
            row["rescoring_features"]["spectralis_score"] = None


def reformat_engine_column(row):
    scan, engines = row["scans"].split("||")

    engine_list = engines.split("|")
    engine_prediction = {}

    for engine in ENGINES:
        if engine in engine_list:
            engine_prediction[engine] = True
        else:
            engine_prediction[engine] = False

    return scan, engine_prediction
