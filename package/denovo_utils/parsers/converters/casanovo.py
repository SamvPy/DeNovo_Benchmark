"""Casanovo parser function."""

import logging
import os

import pandas as pd
from psm_utils import PSM, PSMList
from pyteomics import mgf
from pyteomics.mztab import MzTab
from tqdm import tqdm

from ...utils.proforma import parse_peptidoform

tqdm.pandas()

logger = logging.getLogger(__name__)
logging.basicConfig(filename="denovo_output_parsing.log", level=logging.INFO)


def casanovo_parser(
    result_path: str, mgf_path: str, mapping: dict, max_length: int = 30, **kwargs
) -> PSMList:
    """
    Return a `PSMList` from a casanovo search result and its associated spectral file.

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
    result_path = os.path.splitext(result_path)[0] + ".mztab"

    # ASSUMPTION:
    # The CasaNovo spectra_ref columns has prefix
    # ms_run[1]:index= and the number is a count
    # of spectra, starting from 0 and going to n for a file with n spectra
    mgf_file = pd.DataFrame(pd.DataFrame(mgf.read(mgf_path))["params"].tolist())
    _ = mgf_file.pop("charge")
    result = (
        pd.DataFrame(MzTab(result_path).spectrum_match_table)
        .set_index("PSM_ID")
        .reset_index()
        .reset_index()
    )
    run = os.path.basename(result_path)

    mgf_file = mgf_file.reset_index()
    result["index"] = result.spectra_ref.apply(lambda x: int(x.split("=")[-1]))

    # Fuse the metadata of the spectra with result file
    joined_file = mgf_file.merge(result, on="index")

    # Sanity check
    assert len(result) == len(joined_file)
    if len(mgf_file) > len(joined_file):
        logging.info(
            f"{result_path} - ContraNovo:"
            f" dropped {len(mgf_file)-len(joined_file)} spectra."
        )

    # Parse to psm utils type format
    joined_file["peptidoform"] = joined_file.apply(
        lambda x: str(x["sequence"]) + "/" + str(int(x["charge"])), axis=1
    )
    joined_file["precursor_mz"] = joined_file["pepmass"].apply(lambda x: x[0])

    # Parse peptidoforms with max_length == 30
    joined_file["peptidoform"] = joined_file["peptidoform"].apply(
        lambda x: parse_peptidoform(x, mapping, max_length)
    )
    joined_file = joined_file.dropna(subset=["peptidoform"]).reset_index(drop=True)

    psm_list = joined_file.progress_apply(
        lambda x: PSM(
            peptidoform=x["peptidoform"],
            spectrum_id=x["title"],
            run=run,
            score=x["search_engine_score[1]"],
            precursor_mz=x["precursor_mz"],
            retention_time=x["rtinseconds"],
            source=x["search_engine"][0] + x["search_engine"][1],
            metadata={
                "aa_scores": x["opt_ms_run[1]_aa_scores"],
                "calc_mass_to_charge": x["calc_mass_to_charge"],
                "spectra_ref": x["spectra_ref"],
                "scans": x["scans"],
            },
        ),
        axis=1,
    ).tolist()

    psmlist = PSMList(psm_list=psm_list)
    return psmlist
