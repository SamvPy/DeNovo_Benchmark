"""Spectralis parser class."""

import os

import numpy as np
import pandas as pd
from psm_utils import PSM, PSMList
from psm_utils.io import read_file
from tqdm import tqdm

from ..constants import ENGINES, ENGINES_MAPPING, EXTENSIONS
from .base import DenovoEngineConverter

tqdm.pandas()

## LEGACY

class SpectralisParser:
    """SpectralisParser class."""

    def __init__(self, mgf_path: str, results_dir: str) -> None:
        """
        Spectralis output parser towards PSMList.

        The class structure makes it possible to add results from disparate
        spectralis outputs to a single psmlist, making it easier to analyse.

        Parameters
        ----------
        mgf_path: str
            Path towards the mgf file analysed with spectralis.
        result_dir: str
            Path towards the folder storing all results from
            disparate de novo search engines

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
        self.added_results = {engine: False for engine in ENGINES}
        self.psmlists = {}

    def add_xtandem(self, result_path: str) -> None:
        """
        Load X!Tandem result into the spectralis parser object.

        This result file is parsed with 'percolator' setting in psm-utils.

        Parameter
        ---------
        result_path: str
            Path towards the X!Tandem file, processed with percolator.
        """
        self._check_filename(result_path)
        self.psmlists["xtandem"] = read_file(result_path, filetype="percolator")
        self.added_xtandem = True

    def parse(self, path_spectralis, overwrite=False) -> None:
        """
        Load and parse the spectralis-rescored result files.

        The parsing is performed as follows:
        1. The spectralis file is read and the scores-spectrum-id pairs are extracted.
        2. The original result files for engine, which provided some predictions, are
        loaded by looking into the directory `results_dir/<engine>`.
        3. The scores are stored into the psmlists,
        in the `rescoring_features` dictionary.

        To extract the parsed result files, call properties `psmlist` or `dataframe`.

        Parameters
        ----------
        path_spectralis: str
            Path towards the spectralis file.
        overwrite: bool
            If results already loaded for a given engine,
            whether to overwrite those results or not.
        """
        _ = ["Spectralis_score", "scans", "filename"]

        self._check_filename(path_spectralis)
        filename_spectralis = os.path.basename(path_spectralis).split(".")[0]

        spectralis_out = pd.read_csv(path_spectralis)

        spectralis_out_engines = spectralis_out.apply(reformat_engine_column, axis=1)
        scan_engine_df = pd.DataFrame(
            spectralis_out_engines.tolist(), columns=["scans", "engine"]
        )
        spectralis_out_processed = pd.concat(
            [
                spectralis_out["Spectralis_score"],
                scan_engine_df["scans"],
                pd.DataFrame(scan_engine_df["engine"].tolist()),
            ],
            axis=1,
        )

        spectralis_out_processed, loaded_engines = self._clean_dataframe(
            spectralis_out_processed
        )
        spectralis_out_processed["filename"] = filename_spectralis
        # print(loaded_engines)

        for engine in loaded_engines:
            if not overwrite and self.added_results[engine]:
                print(f"Already loaded results for {engine}! Skipping...")
                continue

            parser = DenovoEngineConverter.select(ENGINES_MAPPING[engine])
            psm_df = parser.parse(
                result_path=os.path.join(
                    self.results_dir,
                    ENGINES_MAPPING[engine],
                    self.filename + EXTENSIONS[ENGINES_MAPPING[engine]],
                ),
                mgf_path=self.mgf_path,
            ).to_dataframe()

            spectralis_engine_results = spectralis_out_processed.loc[
                spectralis_out_processed[engine], ["scans", "Spectralis_score"]
            ].set_index("scans")

            _ = psm_df.progress_apply(
                lambda x: self._add_spectralis_score(x, spectralis_engine_results),
                axis=1,
            )

            self.psmlists[ENGINES_MAPPING[engine]] = PSMList(
                psm_list=[PSM(**x) for x in psm_df.to_dict("records")]
            )
            self.added_results[engine] = True

    @property
    def psmlist(self) -> PSMList:
        """Merge and return the loaded psmlists (parse)."""
        psmlist_all = PSMList(psm_list=[])
        for psmlist in self.psmlists.values():
            psmlist_all += psmlist
        return psmlist_all

    @property
    def dataframe(self) -> pd.DataFrame:
        """Merge and return the loaded psmlists (parse) as dataframe."""
        return self.psmlist.to_dataframe()

    def _check_filename(self, p: str) -> None:
        """
        Check whether the filename matches with the name of the mgf file.

        Parameter
        ---------
        p: str
            Path to some file which will be compared to the filename attribute.
        """
        p_filename = os.path.basename(p).split(".")[0]
        if self.filename not in p_filename:
            print("Warning! Are you sure you loaded the correct result file?")
            print(f"Loaded {p_filename} while mgf is {self.filename}.")

    def _clean_dataframe(self, df: pd.DataFrame):
        """Delete columns (de novo search engine column) which have no entries."""
        engine_empty = np.array(ENGINES)[
            (df.loc[:, ENGINES].sum() == 0).tolist()
        ].tolist()

        loaded_engines = []
        for engine in ENGINES:
            if engine not in engine_empty:
                loaded_engines.append(engine)

        for engine in engine_empty:
            _ = df.pop(engine)
        return df, loaded_engines

    def _add_spectralis_score(self, row, spectralis_score_df):
        """Add spectralis score to the rescoring_features dictionary."""
        score = spectralis_score_df.loc[row["spectrum_id"], "Spectralis_score"]

        # TODO: If I want to include sequences not scored, this should return None
        try:
            row["rescoring_features"]["spectralis_score"] = float(score)
        except KeyError:
            row["rescoring_features"]["spectralis_score"] = None


def reformat_engine_column(row):
    """
    Reformat the engine column.

    The scan column in the spectralis file are parsed as scan|engine1||engine2...
    They are now parsed to a dictionary with booleans.

    Parameter
    ---------
    row: pd.Series
        A series (or dict) with the `scans` column

    Returns
    -------
    scan: str
        The scan id of the spectrum
    engine_prediction: dict
        A boolean dictionary indicating which engine provided the prediction.
    """
    scan, engines = row["scans"].split("||")

    engine_list = engines.split("|")
    engine_prediction = {}

    for engine in ENGINES:
        if engine in engine_list:
            engine_prediction[engine] = True
        else:
            engine_prediction[engine] = False

    return scan, engine_prediction
