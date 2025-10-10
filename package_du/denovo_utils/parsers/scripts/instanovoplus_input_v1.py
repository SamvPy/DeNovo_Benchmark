import argparse
import json
import logging
import os

import pandas as pd
from pyteomics import mgf
from glob import glob

from ..constants import UNSUPPORTED_MODS_INSTANOVO_PLUS, EXTENSIONS
from ..converters import DenovoEngineConverter
from ..exceptions import DenovoEngineNotSupported

logger = logging.getLogger(__name__)
logging.basicConfig(filename="instanovoplust_input_parsing.log", level=logging.INFO)


from denovo_utils.parsers.converters import DenovoEngineConverter
from denovo_utils.parsers.constants import UNSUPPORTED_MODS_INSTANOVO_PLUS
from pyteomics import mgf
import numpy as np
import pandas as pd
import json
import os
import re

class InstanovoPlusParser:
    def __init__(self, engine, mgf_path, result_path):
        self.parser = DenovoEngineConverter.select(engine)
        self.engine = engine
        self.mgf_path = mgf_path
        self.result_path = result_path
        self.tokenizer_regex = (
            # First capture group: matches either:
            # - A UNIMOD annotation like [UNIMOD:123]
            # - Any text inside parentheses like (ox) or (+.98)
            r"(\[UNIMOD:\d+\]|\([^)]+\))|"
            # Second capture group: starts with a valid amino acid letter
            # (including U for selenocysteine and O for pyrrolysine)
            r"([A-Z]"
            # Optionally followed by a UNIMOD annotation
            r"(?:\[UNIMOD:\d+\]|"
            # Or optionally followed by text in parentheses
            r"\([^)]+\))?"
            # Close second capture group
            r")"
        )

    def _parse_input_spectrum(self, psmlist):
        df = psmlist.to_dataframe().reset_index(names='psm_id')
        df["sequence"] = df["peptidoform"].apply(lambda x: x.proforma)
        df["sequence"] = df["sequence"].apply(self._proforma_to_instanovo_mod)
        df["precursor_charge"] = df["peptidoform"].apply(lambda x: x.precursor_charge)

        df = df.rename(columns={"spectrum_id": "title"}).loc[
            :,
            [
                "psm_id",
                "title",
                "sequence",
                "precursor_charge",
                "precursor_mz",
                "rank",
                "source"
            ],
        ]

        mgf_df = pd.DataFrame(mgf.read(self.mgf_path))
        mgf_df = pd.concat(
            [
                pd.DataFrame(mgf_df["params"].tolist()),
                mgf_df[["m/z array", "intensity array"]],
            ],
            axis=1,
        )[["title", "m/z array", "intensity array"]]

        df = df.merge(mgf_df, on="title")
        df['specid'] = df['title']
        df["title"] = df["title"].apply(lambda x: x+f"||{self.engine}")

        df = df.reset_index()
        return df

    def _parse_input_predictions(self, df):
        input_preds = df[['psm_id', 'specid', 'sequence']].copy()
        input_preds['predictions'] = input_preds['sequence']
        input_preds['predictions_tokenised'] = input_preds['predictions'].apply(self._tokenize)
        input_preds['log_probabilities'] = np.nan
        input_preds['token_log_probs_col'] = np.nan
        return input_preds


    def parse(self):
        psmlist = self.parser.parse(
            result_path=self.result_path,
            mgf_path=self.mgf_path
        )
        if len(psmlist) == 0:
            return False

        df = self._parse_input_spectrum(psmlist)
        n_before_filter = len(df)
        df = df.dropna(subset='sequence')
        n_after_filter = len(df)
        if n_before_filter != n_after_filter:
            logging.warning(
                "Filtered {} rows due to empty peptide predictions!".format(
                    n_before_filter-n_after_filter
                )
            )
        self.index_map = df[["index", "title"]].set_index("index").to_dict("index")
        
        df = df.loc[
            :,
            [
                "psm_id",
                "title",
                "specid",
                "rank",
                "index",
                "sequence",
                "precursor_mz",
                "precursor_charge",
                "m/z array",
                "intensity array",
                "source"
            ],
        ].rename(columns={
            "m/z array": "mz_array",
            "intensity array": "intensity_array",
            'index': 'scan_number',
        })
        self.df = df.reset_index(drop=True)

        self.intitial_predictions = self._parse_input_predictions(
            df
        )
        return True

    
    def write(self, filename = 'filename', out_dir = '.'):
        self.df.to_feather(os.path.join(out_dir, f'{filename}.ipc'))
        self.intitial_predictions.to_csv(os.path.join(out_dir, f'{filename}.csv'))

    def store_dictionary_to_json(self, dictionary: dict, filename: str):
        with open(filename, "w") as json_file:
            json.dump(dictionary, json_file, indent=4)

    def _proforma_to_instanovo_mod(self, peptidoform):
        for not_supported_mod_label in UNSUPPORTED_MODS_INSTANOVO_PLUS:
            if not_supported_mod_label in peptidoform:
                return None
        return peptidoform.split("/")[0]
    
    def _tokenize(self, peptide):
        if peptide is None:
            return None
        return [
            item
            for sublist in re.findall(self.tokenizer_regex, peptide)
            for item in sublist
            if item
        ]


def main(args):
    in_parser = InstanovoPlusParser(
        engine=args.denovo_engine,
        mgf_path=args.mgf_path,
        result_path=args.result_file
    )
    parsed = in_parser.parse()
    if parsed:
        in_parser.write(os.path.basename(args.result_file).split('.')[0])

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Write mgf files, annotated by a given de novo search"
    )

    parser.add_argument(
        "-m", "--mgf_path", required=True, help="Path to unannotated mgf-file"
    )
    parser.add_argument(
        "-d",
        "--denovo_engine",
        required=True,
        help="The denovo engine used to generate the files search result file.",
    )
    parser.add_argument(
        "-r",
        "--result_file",
        required=True,
        help="Result file from a database or de novo search engine to refine."
    )
    parser.add_argument(
        "-o",
        "--output_folder",
        default="",
        help="Output folder to store annotated mgf-file",
    )
    parser.add_argument(
        "-x",
        "--exclusion_list",
        default=[],
        help="List containing tags in peptides that should be dropped (e. g. modification tags...)",
    )

    args = parser.parse_args()

    main(args)
