import pandas as pd
import os
import logging

from pyteomics import mgf
from psm_utils import PSM, PSMList

from tqdm import tqdm

from .utils import parse_peptidoform

tqdm.pandas()

def novob_parser(result_path: str, mgf_path: str, mapping: dict, max_length=30):
    result_path = os.path.splitext(result_path)[0] + ".tsv"

    mgf_file = pd.DataFrame(pd.DataFrame(mgf.read(mgf_path))["params"].tolist())
    _ = mgf_file.pop("charge")
    result = pd.read_csv(result_path, sep="\t", header=None).rename(
        columns={
            0: "Mcount",
            1: "charge",
            2: "peptide_mass",
            3: "sequence_forward",
            4: "mass_forward",
            5: "probability_forward",
            6: "sequence_reverse",
            7: "mass_reverse",
            8: "probability_reverse",
            9: "scans"
        }
    )
    result["scans"] = result.apply(lambda x: x["scans"][2:-1], axis=1) # Remove b' ... '
    run = os.path.basename(result_path)

    # Fuse the metadata of the spectra with result file
    joined_file = result.merge(mgf_file, on="scans")

    # Sanity checks
    assert len(result) == len(joined_file)
    if len(mgf_file) > len(joined_file):
        logging.info(f"{result_path} for NovoB: dropped {len(mgf_file)-len(joined_file)} spectra.")

    # Parse to psm utils type format
    joined_file["precursor_mz"] = joined_file["pepmass"].apply(
        lambda x: x[0]
    )

    joined_file = joined_file.progress_apply(
        lambda x: select_top_PSM(x, max_length, mapping, run),
        axis=1
    )
    joined_file = joined_file.dropna()

    psmlist = PSMList(
        psm_list=joined_file.tolist()
    )
    return psmlist


def select_top_PSM(x, max_length, mapping, run=None):
    best = ""

    if x["probability_forward"] > x["probability_reverse"]:
        best = "forward"
    elif x["probability_forward"] < x["probability_reverse"]:
        best = "reverse"
    elif abs(x["mass_forward"]) < abs(x["mass_reverse"]):
        best = "forward"
    else:
        best = "reverse"

    peptidoform = eval(x[f"sequence_{best}"])[0] + "/" + str(int(x["charge"]))
    mass_error = x[f"mass_{best}"]
    proba = x[f"probability_{best}"]

    peptide = parse_peptidoform(peptidoform, mapping, max_length)
    if isinstance(peptide, type(None)):
        return None
    else:
        return PSM(
            peptidoform=peptide,
            spectrum_id=x["title"],
            run=run,
            score=proba,
            precursor_mz=x["precursor_mz"],
            retention_time=x["rtinseconds"],
            source="NovoB",
            metadata={
                "ppm_error": mass_error,
                "scans": x["scans"]
            }
        )