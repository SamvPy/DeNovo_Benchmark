"""DeepNovo parser function."""

import os

import pandas as pd
from psm_utils import PSM, PSMList
from pyteomics import mgf
from tqdm import tqdm

from ...utils.proforma import parse_peptidoform
from .utils import mzml_reader

tqdm.pandas()


def deepnovo_parser(
    result_path: str, mgf_path: str, mapping: dict, max_length=30, **kwargs
) -> PSMList:
    """
    Return a `PSMList` from a DeepNovo search result and its MGF spectral file.

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

    if mgf_path.lower().endswith('.mzml'):
        mgf_file = mzml_reader(mgf_path)
    else:
        mgf_file = pd.DataFrame(pd.DataFrame(mgf.read(mgf_path))["params"].tolist())

    result = pd.read_csv(result_path, sep='\t')[['scan', "output_seq", "output_score"]]
    result = result.rename(columns={'scan': 'scans'})
    result['scans'] = result['scans'].apply(str)
    run = os.path.basename(result_path)

    # Fuse the metadata of the spectra with result file
    joined_file = result.merge(mgf_file, on='scans')

    # Sanity checks
    # assert len(result) == len(joined_file)
    joined_file = joined_file.dropna(subset=["output_seq"]).reset_index(drop=True)

    joined_file['peptidoform'] = joined_file.apply(
        lambda x: "".join(x['output_seq'].split(',')) + "/" + str(int(x["charge"][0])),
        axis=1
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
                score=x["output_score"],
                precursor_mz=x["precursor_mz"],
                retention_time=x["rtinseconds"]/60,
                source="DeepNovo",
                metadata={},
            ),
            axis=1,
        ).tolist()
    )
    return psmlist