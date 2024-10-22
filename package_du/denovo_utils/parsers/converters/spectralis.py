"""Spectralis parser class."""

import os

import pandas as pd
from pyteomics import mgf
from psm_utils import PSM, PSMList
from tqdm import tqdm

from ...utils.proforma import parse_peptidoform

tqdm.pandas()


def spectralis_parser(
    result_path: str, mgf_path: str, mapping: dict, max_length: int = 30, **kwargs
) -> PSMList:
    
    result_path = os.path.splitext(result_path)[0] + ".csv"
    
    mgf_file = pd.DataFrame(pd.DataFrame(mgf.read(mgf_path))["params"].tolist())
    result = pd.read_csv(result_path)
    run = os.path.basename(result_path)


    result['source'] = result['scans'].apply(lambda x: x.split("||")[1])
    result['title'] = result['scans'].apply(lambda x: x.split("||")[0])

    joined_file = result.merge(mgf_file, on="title")

    # Sanity check
    assert len(joined_file) == len(result)


    joined_file["peptidoform"] = joined_file.apply(
        lambda x: str(x["peptide_spectralis-ea"]) + "/" + str(int((x["charge"][0]))), axis=1
    )
    joined_file["peptidoform"] = joined_file["peptidoform"].apply(
        lambda x: parse_peptidoform(x, mapping, max_length)
    )
    joined_file["peptidoform_init"] = joined_file.apply(
        lambda x: str(x["peptide_init"]) + "/" + str(int((x["charge"][0]))), axis=1
    )
    joined_file["peptidoform_init"] = joined_file["peptidoform_init"].apply(
        lambda x: parse_peptidoform(x, mapping, max_length)
    )

    joined_file["precursor_mz"] = joined_file["pepmass"].apply(lambda x: x[0])

    psmlist = PSMList(
        psm_list=joined_file.progress_apply(
            lambda x: PSM(
                peptidoform=x["peptidoform"],
                spectrum_id=x["title"],
                run=run,
                score=x["score_spectralis-ea"],
                precursor_mz=x["precursor_mz"],
                retention_time=x["rtinseconds"],
                source="Spectralis",
                metadata={
                    "base_prediction": x["source"]
                },
                provenance_data={
                    "init_score_spectralis": x["score_init"],
                    "init_peptidoform": x["peptidoform_init"].proforma}
            ),
            axis=1,
        ).tolist()
    )
    return psmlist
