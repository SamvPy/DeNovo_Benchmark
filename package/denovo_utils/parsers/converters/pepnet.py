import pandas as pd
import os
import logging

from pyteomics.mztab import MzTab
from pyteomics import mgf
from psm_utils import PSM, PSMList

from tqdm import tqdm

from .utils import parse_peptidoform

tqdm.pandas()

def pepnet_parser(result_path: str, mgf_path: str, mapping: dict, max_length=30):
    result_path = os.path.splitext(result_path)[0] + ".tsv"

    mgf_file = pd.DataFrame(pd.DataFrame(mgf.read(mgf_path))["params"].tolist())
    result = pd.read_csv(result_path, sep="\t").rename(
        columns={"TITLE": "title"}
    )
    run = os.path.basename(result_path)

    # Fuse the metadata of the spectra with result file
    joined_file = result.merge(mgf_file, on="title")

    # Sanity checks
    assert len(result) == len(joined_file)

    joined_file["peptidoform"] = joined_file.apply(
        lambda x: str(x["DENOVO"]) + "/" + str(int((x["charge"][0]))), axis=1
    )
    joined_file["precursor_mz"] = joined_file["pepmass"].apply(
        lambda x: x[0]
    )

    joined_file["peptidoform"] = joined_file["peptidoform"].apply(
        lambda x: parse_peptidoform(x, mapping, max_length)
    )
    joined_file = joined_file.dropna(subset=["peptidoform"]).reset_index(drop=True)

    psmlist = PSMList(
        psm_list=joined_file.progress_apply(
            lambda x: PSM(
                peptidoform=x["peptidoform"],
                spectrum_id=x["title"],
                run=run,
                score=x["Score"],
                precursor_mz=x["precursor_mz"],
                retention_time=x["rtinseconds"],
                source="PepNet",
                metadata={
                    "positional_scores": x["Positional Score"],
                    "ppm_error": x["PPM Difference"]
                }
            ),
            axis=1
        ).tolist()
    )
    return psmlist