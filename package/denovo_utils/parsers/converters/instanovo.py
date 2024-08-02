import pandas as pd
import os

from pyteomics import mgf
from psm_utils import PSM, PSMList

from tqdm import tqdm

from ...utils.proforma import parse_peptidoform

tqdm.pandas()

def instanovo_parser(result_path: str, mgf_path: str, mapping: dict, max_length=30, **kwargs):
    result_path = os.path.splitext(result_path)[0] + ".csv"


    # ASSUMPTION: 
    # The output of Instanovo has the same length and order (in terms of spectra) as the mgf-file
    mgf_file = pd.DataFrame(pd.DataFrame(mgf.read(mgf_path))["params"].tolist())
    result = pd.read_csv(result_path)
    run = os.path.basename(result_path)

    length_mgf = len(mgf_file)
    length_result = len(result)
    assert length_mgf == length_result

    # Fuse the metadata of the spectra with result file
    joined_file = pd.concat([mgf_file, result], axis=1)
    length_join = len(joined_file)

    # Sanity checks
    assert length_mgf == length_join
    assert length_result == length_join

    joined_file["peptidoform"] = joined_file.apply(
        lambda x: str(x["preds"]) + "/" + str(int((x["charge"][0]))), axis=1
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
                score=x["log_probs"],
                precursor_mz=x["precursor_mz"],
                retention_time=x["rtinseconds"],
                source="InstaNovo",
                metadata={
                    "scans": x["scans"]
                }
            ),
            axis=1
        ).tolist()
    )
    return psmlist