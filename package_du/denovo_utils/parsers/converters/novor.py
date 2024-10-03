"""Novor parser function."""

import logging
import os

import pandas as pd
from psm_utils import PSM, PSMList
from pyteomics import mgf
from tqdm import tqdm

from ...utils.proforma import parse_peptidoform

tqdm.pandas()

logger = logging.getLogger(__name__)
logging.basicConfig(filename="denovo_output_parsing.log", level=logging.INFO)


def novor_parser(
    result_path: str, mgf_path: str, mapping: dict, max_length=30, **kwargs
) -> PSMList:
    """
    Return a `PSMList` from a ContraNovo search result and its MGF spectral file.

    Parameters
    ----------
    result_path: str
        Path to the search results.
    mgf_path: str
        Path to the MGF-file that was used to generate search results
    mapping: dict
        Mapping dictionary for converting engine-specific to
        unified modification and amino acid labels.
    max_length: int
        Any peptide sequence longer than this value will be ignored.

    Return
    ------
    psm_utils.PSMList
        The PSMList, representing the search results.
    """
    result_path = os.path.join(
        os.path.dirname(result_path),
        os.path.basename(result_path).split(".")[0] + ".novor.csv",
    )
    run = os.path.basename(result_path)

    mgf_file = pd.DataFrame(
        pd.DataFrame(mgf.read(mgf_path))["params"].tolist()
    ).reset_index()
    result = pd.read_csv(result_path, sep=", ", skiprows=20, engine="python").rename(
        columns={"scanNum": "scans", "aaScore,": "aaScore"}
    )
    result["scans"] = result["scans"].apply(lambda x: str(x))
    result["index"] = result["# id"].apply(lambda x: int(x) - 1)

    logging.info(
        f"Diff: ({mgf_file.shape[0]-result.shape[0]});"
        f"\tResult shape ({result.shape}), mgf shape ({mgf_file.shape})"
    )

    # joined_file = result.merge(mgf_file, on="scans")
    joined_file = result.merge(mgf_file, on="index")

    assert len(result) == len(joined_file)
    if len(mgf_file) > len(joined_file):
        logging.info(
            f"{result_path} for Novor:"
            f" dropped {len(mgf_file)-len(joined_file)} spectra."
        )

    joined_file["peptidoform"] = joined_file.apply(
        lambda x: str(x["peptide"] + "/" + str(x["z"])), axis=1
    )
    joined_file["peptidoform"] = joined_file["peptidoform"].apply(
        lambda x: parse_peptidoform(x, mapping=mapping, max_length=max_length)
    )
    joined_file["precursor_mz"] = joined_file["pepmass"].apply(lambda x: x[0])

    joined_file = joined_file.loc[
        :,
        [
            "peptidoform",
            "title",
            "score",
            "precursor_mz",
            "rtinseconds",
            "scans_x",
            "err(data-denovo)",
            "ppm(1e6*err/(mz*z))",
            "aaScore",
        ],
    ]
    joined_file = joined_file.dropna(subset=["peptidoform"]).reset_index(drop=True)

    psm_list = joined_file.progress_apply(
        lambda x: PSM(
            peptidoform=x["peptidoform"],
            spectrum_id=x["title"],
            run=run,
            score=x["score"],
            precursor_mz=x["precursor_mz"],
            retention_time=x["rtinseconds"],
            source="Novor",
            metadata={
                "scans": x["scans_x"],
                "err(data-denovo)": x["err(data-denovo)"],
                "ppm_error": x["ppm(1e6*err/(mz*z))"],
                "aaScore": x["aaScore"],
            },
        ),
        axis=1,
    ).tolist()

    psmlist = PSMList(psm_list=psm_list)
    return psmlist
