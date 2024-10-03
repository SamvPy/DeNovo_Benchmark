import argparse
import json
import logging
import os

import pandas as pd
from pyteomics import mgf
from glob import glob

from ..constants import UNSUPPORTED_MODS_INSTANOVO_PLUS
from ..converters import DenovoEngineConverter
from ..exceptions import DenovoEngineNotSupported

logger = logging.getLogger(__name__)
logging.basicConfig(filename="instanovoplust_input_parsing.log", level=logging.INFO)


def proforma_to_instanovo_mod(peptidoform: str):
    proforma_to_instanovo = {
        "M[UNIMOD:35]": "M(+15.99)",
        "N[UNIMOD:7]": "N(+.98)",
        "Q[UNIMOD:7]": "Q(+.98)",
        "C[UNIMOD:4]": "C(+57.02)",
    }

    for proforma_mod, instanovo_mod in proforma_to_instanovo.items():
        peptidoform = peptidoform.replace(proforma_mod, instanovo_mod)


    for not_supported_mod_label in UNSUPPORTED_MODS_INSTANOVO_PLUS:
        if not_supported_mod_label in peptidoform:
            return None
    return peptidoform.split("/")[0]


def store_dictionary_to_json(dictionary: dict, filename: str):
    with open(filename, "w") as json_file:
        json.dump(dictionary, json_file, indent=4)


def main(args):

    filename = os.path.basename(args.mgf_path).split(".")[0]
    df_list = []

    for filepath in glob(f"./{filename}.*"):
        if filepath.endswith(".mgf"):
            continue
        engine = os.path.basename(filepath).split(".")[1]

        try:
            denovo_engine = DenovoEngineConverter.select(label=engine)
        except DenovoEngineNotSupported as err:
            logging.info(err)
            logging.info("Skipping {}...".format(filepath))
            continue

        psmlist = denovo_engine.parse(result_path=filepath, mgf_path=args.mgf_path)

        df = psmlist.to_dataframe()
        df["sequence"] = df["peptidoform"].apply(lambda x: x.sequence)
        df["modified_sequence"] = df["peptidoform"].apply(lambda x: x.proforma)
        df["modified_sequence"] = df["modified_sequence"].apply(proforma_to_instanovo_mod)
        df["precursor_charge"] = df["peptidoform"].apply(lambda x: x.precursor_charge)

        df = df.rename(columns={"spectrum_id": "title"}).loc[
            :,
            [
                "title",
                "sequence",
                "modified_sequence",
                "precursor_mz",
                "precursor_charge",
            ],
        ]

        mgf_df = pd.DataFrame(mgf.read(args.mgf_path))
        mgf_df = pd.concat(
            [
                pd.DataFrame(mgf_df["params"].tolist()),
                mgf_df[["m/z array", "intensity array"]],
            ],
            axis=1,
        )[["title", "m/z array", "intensity array"]]

        df = df.merge(mgf_df, on="title")
        df["title"] = df["title"].apply(lambda x: x+f"||{engine}")

        df_list.append(df)

    dfs = pd.concat(df_list, ignore_index=True)

    dfs = dfs.reset_index()
    index_map = dfs[["index", "title"]].set_index("index").to_dict("index")
    store_dictionary_to_json(index_map, filename + ".json")

    dfs = dfs.loc[
        :,
        [
            "index",
            "sequence",
            "modified_sequence",
            "precursor_mz",
            "precursor_charge",
            "m/z array",
            "intensity array",
        ],
    ].rename(columns={"m/z array": "mz_array", "intensity array": "intensity_array"})


    dfs.to_feather(os.path.join(args.output_folder, filename) + ".feather")


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
