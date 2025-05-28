"""InstaNovo parser function."""

import os

import pandas as pd
from psm_utils import PSM, PSMList
from pyteomics import mgf
from tqdm import tqdm

from ...utils.proforma import parse_peptidoform
from .utils import mzml_reader, infer_rank

tqdm.pandas()


def instanovo_parser(
    result_path: str, mgf_path: str, mapping: dict, max_length=30, **kwargs
):
    """
    Return a `PSMList` from a InstaNovo search result and its MGF spectral file.

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
    result_path = os.path.splitext(result_path)[0] + ".csv"

    # ASSUMPTION:
    # The output of Instanovo has the same length and order
    # (in terms of spectra) as the mgf-file
    if mgf_path.lower().endswith('.mzml'):
        mgf_file = mzml_reader(mgf_path)
    else:
        mgf_file = pd.DataFrame(pd.DataFrame(mgf.read(mgf_path))["params"].tolist())

    result = pd.read_csv(result_path)
    _ = result.pop('precursor_mz')
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
        lambda x: str(x["predictions"]) + "/" + str(int((x["precursor_charge"]))), axis=1
    )
    joined_file["precursor_mz"] = joined_file["pepmass"].apply(lambda x: x[0])
    joined_file["peptidoform"] = joined_file["peptidoform"].apply(
        lambda x: parse_peptidoform(x, mapping, max_length)
    )
    joined_file = joined_file.dropna(subset=["peptidoform"]).reset_index(drop=True)

    joined_file['rank'] = joined_file.groupby('title')['log_probabilities'].rank(
        ascending=False, method='dense'
    )

    tqdm.pandas(desc='Parsing Instanovo results to PSMList')
    psmlist = PSMList(
        psm_list=joined_file.progress_apply(
            lambda x: PSM(
                peptidoform=x["peptidoform"],
                spectrum_id=x["title"],
                run=run,
                rank=x['rank'],
                score=x["log_probabilities"],
                precursor_mz=x["precursor_mz"],
                retention_time=x["rtinseconds"]/60,
                source="InstaNovo",
                metadata={
                    "scan_number": x["scan_number"],
                    "aa_scores": x["token_log_probabilities"]
                },
            ),
            axis=1,
        ).tolist()
    )
    return psmlist

def instanovoplus_parser(
    result_path: str, mgf_path: str, mapping: dict, max_length=30, **kwargs
):
    result_path = os.path.splitext(result_path)[0] + ".csv"

    # ASSUMPTION:
    # The output of Instanovo has the same length and order
    # (in terms of spectra) as the mgf-file
    if mgf_path.lower().endswith('.mzml'):
        mgf_file = mzml_reader(mgf_path)
    else:
        mgf_file = pd.DataFrame(pd.DataFrame(mgf.read(mgf_path))["params"].tolist())

    result = pd.read_csv(result_path)
    run = os.path.basename(result_path)

    # Fuse the metadata of the spectra with result file
    result["title"] = result["specid"]
    joined_file = result.merge(mgf_file, on="title")

    # Sanity check
    assert len(joined_file) == len(result)

    joined_file["peptidoform"] = joined_file.apply(
        lambda x: str(x["diffusion_predictions"]) + "/" + str(int(x["precursor_charge"])), axis=1
    )
    joined_file["precursor_mz"] = joined_file["pepmass"].apply(lambda x: x[0])
    joined_file["peptidoform"] = joined_file["peptidoform"].apply(
        lambda x: parse_peptidoform(x, mapping, max_length)
    )
    joined_file = joined_file.dropna(subset=["peptidoform"]).reset_index(drop=True)

    tqdm.pandas(desc='Parsing Instanovo+ results to PSMList')
    psmlist = PSMList(
        psm_list=joined_file.progress_apply(
            lambda x: PSM(
                peptidoform=x["peptidoform"],
                spectrum_id=x["specid"],
                run=run,
                rank=x['rank'],
                score=x["diffusion_log_probabilities"],
                precursor_mz=x["precursor_mz"],
                retention_time=x["rtinseconds"]/60,
                source="InstaNovo+",
                metadata={
                    "base_prediction": x["source"],
                    "base_peptidoform":str(x["transformer_predictions"]) + "/" + str(int(x["precursor_charge"]))},
            ),
            axis=1,
        ).tolist()
    )
    return psmlist


def instanovo_parser_legacy(
    result_path: str, mgf_path: str, mapping: dict, max_length=30, **kwargs
):
    """
    Return a `PSMList` from a InstaNovo search result and its MGF spectral file.

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
    result_path = os.path.splitext(result_path)[0] + ".csv"

    # ASSUMPTION:
    # The output of Instanovo has the same length and order
    # (in terms of spectra) as the mgf-file
    if mgf_path.lower().endswith('.mzml'):
        mgf_file = mzml_reader(mgf_path)
    else:
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
    joined_file["precursor_mz"] = joined_file["pepmass"].apply(lambda x: x[0])
    joined_file["peptidoform"] = joined_file["peptidoform"].apply(
        lambda x: parse_peptidoform(x, mapping, max_length)
    )
    joined_file = joined_file.dropna(subset=["peptidoform"]).reset_index(drop=True)

    joined_file['rank'] = joined_file.groupby('title')['log_probs'].rank(
        ascending=False, method='dense'
    )

    tqdm.pandas(desc='Parsing Instanovo results to PSMList')
    psmlist = PSMList(
        psm_list=joined_file.progress_apply(
            lambda x: PSM(
                peptidoform=x["peptidoform"],
                spectrum_id=x["title"],
                run=run,
                rank=x['rank'],
                score=x["log_probs"],
                precursor_mz=x["precursor_mz"],
                retention_time=x["rtinseconds"]/60,
                source="InstaNovo",
                metadata={"scans": x["scans"]},
            ),
            axis=1,
        ).tolist()
    )
    return psmlist

def instanovoplus_parser_legacy(
    result_path: str, mgf_path: str, mapping: dict, max_length=30, **kwargs
):
    result_path = os.path.splitext(result_path)[0] + ".csv"

    # ASSUMPTION:
    # The output of Instanovo has the same length and order
    # (in terms of spectra) as the mgf-file
    if mgf_path.lower().endswith('.mzml'):
        mgf_file = mzml_reader(mgf_path)
    else:
        mgf_file = pd.DataFrame(pd.DataFrame(mgf.read(mgf_path))["params"].tolist())

    result = pd.read_csv(result_path)
    run = os.path.basename(result_path)

    # Fuse the metadata of the spectra with result file
    result["title"] = result["spectrum_id"].apply(lambda x: x.split("||")[0])
    result["source"] = result["spectrum_id"].apply(lambda x: x.split("||")[1])
    joined_file = result.merge(mgf_file, on="title")

    # Sanity check
    assert len(joined_file) == len(result)

    joined_file["peptidoform"] = joined_file.apply(
        lambda x: str(x["predictions"]) + "/" + str(int((x["charge"][0]))), axis=1
    )
    joined_file["precursor_mz"] = joined_file["pepmass"].apply(lambda x: x[0])
    joined_file["peptidoform"] = joined_file["peptidoform"].apply(
        lambda x: parse_peptidoform(x, mapping, max_length)
    )
    joined_file = joined_file.dropna(subset=["peptidoform"]).reset_index(drop=True)

    spectrum_id_dict = {}

    tqdm.pandas(desc='Parsing Instanovo+ results to PSMList')
    psmlist = PSMList(
        psm_list=joined_file.progress_apply(
            lambda x: PSM(
                peptidoform=x["peptidoform"],
                spectrum_id=x["title"],
                run=run,
                rank=infer_rank(x['title'], spectrum_id_dict),
                score=x["log_probabilities"],
                precursor_mz=x["precursor_mz"],
                retention_time=x["rtinseconds"]/60,
                source="InstaNovo+",
                metadata={"base_prediction": x["source"]},
            ),
            axis=1,
        ).tolist()
    )
    return psmlist

