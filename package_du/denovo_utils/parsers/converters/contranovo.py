"""ContraNovo parser function."""

import logging
import os

import pandas as pd
from psm_utils import PSM, PSMList
from pyteomics import mgf
from pyteomics.mztab import MzTab
from tqdm import tqdm

from ...utils.proforma import parse_peptidoform
from .utils import mzml_reader

tqdm.pandas()

logger = logging.getLogger(__name__)
logging.basicConfig(filename="denovo_output_parsing.log", level=logging.INFO)


def contranovo_parser(
    result_path: str, mgf_path: str, mapping: dict, max_length=30, im: bool=False, **kwargs
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
    result_path = os.path.splitext(result_path)[0] + ".mztab"

    if mgf_path.lower().endswith('.mzml'):
        mgf_file = mzml_reader(mgf_path, im=im)
    else:
        mgf_file = pd.DataFrame(pd.DataFrame(mgf.read(mgf_path))["params"].tolist())

    _ = mgf_file.pop("charge")

    try:
        result = (
            pd.DataFrame(MzTab(result_path).spectrum_match_table)
            .set_index("PSM_ID")
            .reset_index()
            .reset_index()
        )
    except Exception:
        return PSMList(psm_list=[])
    run = os.path.basename(result_path)

    # Fuse the metadata of the spectra with result file
    joined_file = result.rename(columns={"spectra_ref": "title"}).merge(
        mgf_file, on="title"
    )

    # Sanity checks
    assert len(result) == len(joined_file)
    if len(mgf_file) > len(joined_file):
        logging.info(
            f"{result_path} for ContraNovo:"
            f" dropped {len(mgf_file)-len(joined_file)} spectra."
        )

    # Parse to psm utils type format
    joined_file["peptidoform"] = joined_file.apply(
        lambda x: str(x["sequence"]) + "/" + str(int((x["charge"]))), axis=1
    )
    joined_file["precursor_mz"] = joined_file["pepmass"].apply(lambda x: x[0])
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
                score=x["search_engine_score[1]"],
                precursor_mz=x["precursor_mz"],
                retention_time=x["rtinseconds"]/60,
                source="ContraNovo",
                metadata={
                    "aa_scores": x["opt_ms_run[1]_aa_scores"],
                    "calc_mass_to_charge": x["calc_mass_to_charge"],
                    "scans": x["scans"],
                },
            ),
            axis=1,
        ).tolist()
    )
    return psmlist
