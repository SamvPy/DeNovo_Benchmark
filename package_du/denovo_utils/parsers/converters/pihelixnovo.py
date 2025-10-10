"""pihelixnovo parser function."""

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


def pihelixnovo_parser(
    result_path: str, mgf_path: str, mapping: dict, max_length: int = 30, im: bool=False, **kwargs
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
    result_path = os.path.splitext(result_path)[0] + ".tsv"

    # ASSUMPTION:
    # The CasaNovo spectra_ref columns has prefix
    # ms_run[1]:index= and the number is a count
    # of spectra, starting from 0 and going to n for a file with n spectra
    # https://github.com/Noble-Lab/casanovo/issues/309
    if mgf_path.lower().endswith('.mzml'):
        mgf_file = mzml_reader(mgf_path, im=im)
        # In mzml format, the indexation with spectra_ref is done with the id (title) part

    else:
        mgf_file = pd.DataFrame(pd.DataFrame(mgf.read(mgf_path))["params"].tolist())

    result = pd.read_csv(result_path, header=None, sep='\t')
    result = result.rename(
        columns={
            0: 'title',
            1: 'prediction',
            2: 'score'
        }
    )
    run = os.path.basename(result_path)

    mgf_file = mgf_file.reset_index()

    # Fuse the metadata of the spectra with result file
    joined_file = mgf_file.merge(result, on="title")

    # Sanity check
    # assert len(result) == len(joined_file)

    if len(mgf_file) > len(joined_file):
        logging.info(
            f"{result_path} - pi-HelixNovo:"
            f" dropped {len(mgf_file)-len(joined_file)} spectra."
        )

    # Parse to psm utils type format
    joined_file['peptidoform'] = joined_file.apply(
        lambda x: parse_peptidoform(
            peptide=f"{x['prediction']}/{str(int((x['charge'][0])))}",
            mapping=mapping,
            max_length=max_length
        ),
        axis=1
    )
    joined_file["precursor_mz"] = joined_file["pepmass"].apply(lambda x: x[0])

    joined_file = joined_file.dropna(subset=["peptidoform"]).reset_index(drop=True)
    if len(joined_file) == 0:
        return PSMList(psm_list=[])

    # Check if ion mobility column is present (only when reading mzml files)
    if "ion_mobility" not in joined_file.columns:
        joined_file['ion_mobility'] = None

    joined_file['rank'] = joined_file.groupby('title')['score'].rank(
        ascending=False, method='dense'
    )

    tqdm.pandas(desc='Parsing pi-HelixNovo results to PSMList')
    psm_list = joined_file.progress_apply(
        lambda x: PSM(
            peptidoform=x["peptidoform"],
            spectrum_id=x["title"],
            run=run,
            rank=x['rank'],
            score=x["score"],
            precursor_mz=x["precursor_mz"],
            retention_time=x["rtinseconds"]/60,
            ion_mobility=x['ion_mobility'],
            source='pi-HelixNovo',
            metadata={},
        ),
        axis=1,
    ).tolist()

    psmlist = PSMList(psm_list=psm_list)
    return psmlist
