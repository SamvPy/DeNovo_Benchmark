from ..constants import PEPNOVO_COLUMN_MAPPING, PEPNOVO_COLUMNS
import re
import pandas as pd
import os
from io import TextIOWrapper
from ...utils.proforma import parse_peptidoform
from pyteomics import mgf
from psm_utils import PSMList, PSM
from tqdm import tqdm
tqdm.pandas()

def parse_metadata_head(input_string) -> tuple:
    # Define the regex pattern
    pattern = r'>>\s*(\d+)\s+(\d+)\s+(.+?)(?:\s+\((SQS\s+\d+\.\d+)\))?$'
    
    # Use re.match to match the pattern
    match = re.match(pattern, input_string)
    
    if match:
        # Extract the matched groups
        number1 = match.group(1)
        number2 = match.group(2)
        random_sequence = match.group(3)
        sqs = match.group(4)
        if not sqs:
            sqs = ""
        return [number1, number2, random_sequence, sqs]
    else:
        raise ValueError(f"The input string does not match the expected format: {input_string}")
    
def parse_header(header):
    header_reformatted = []
    header_list = header.split("\t")
    for col in header_list:
        header_reformatted.append(PEPNOVO_COLUMN_MAPPING[col])
    return header_reformatted

def parse_entry(entry_list: list):
    # Check if a prediction was outputted
    # for error_string, error_type in spectrum_errors.items():
    #     if error_string in entry_list:
    #         return []
    if not entry_list[1].startswith("#Index"):
        return []
    
    metadata_head = entry_list[0]
    # header = entry_list[1]
    inferences = entry_list[2:]
    
    metadata = parse_metadata_head(metadata_head)
    # column_names = parse_header(header)

    # Only get the top prediction
    top_candidate = inferences[0].split("\t")
    entry = top_candidate + metadata
    return entry[1:]


def pepnovo_to_df(file: TextIOWrapper) -> pd.DataFrame:

    entries = []
    entry_list = []
    entry = []

    for line in file:
        line=line.strip()
        if line.startswith(">>"):

            if entry_list:
                entry = parse_entry(entry_list)
            if entry:
                entries.append(entry)
            
            entry_list = [line]
            continue

        entry_list.append(line)

    entry = parse_entry(entry_list)
    if entry:
        entries.append(entry)
    return pd.DataFrame(entries, columns=PEPNOVO_COLUMNS[1:])



def pepnovo_parser(result_path: str, mgf_path: str, mapping: dict, max_length=30, **kwargs):
    result_path = os.path.join(
        os.path.dirname(result_path), os.path.basename(result_path).split(".")[0]+'.mgf.out'
    )
    run = os.path.basename(result_path)
    mgf_file = pd.DataFrame(pd.DataFrame(mgf.read(mgf_path))["params"].tolist())
    # Do not trust MS1-inferred charge state
    # This not matching could be used as a feature?
    _ = mgf_file.pop("charge")

    # Parse the pepnovo text file towards a dataframe
    with open(result_path) as f:
        result = pepnovo_to_df(f)

    # Handle the dataframe as all other parsers
    joined_file = result.merge(mgf_file, on="title")
    assert len(result)==len(joined_file)

    joined_file["peptidoform"] = joined_file.apply(
        lambda x: x["peptide"] + "/" + str(int(x["charge"])), axis=1
    )
    joined_file["precursor_mz"] = joined_file["pepmass"].apply(
        lambda x: x[0]
    )
    joined_file["peptidoform"] = joined_file["peptidoform"].apply(
        lambda x: parse_peptidoform(x, mapping=mapping, max_length=max_length)
    )
    joined_file = joined_file.dropna(subset=["peptidoform"]).reset_index(drop=True)

    psm_list = joined_file.progress_apply(
            lambda x: PSM(
                peptidoform=x["peptidoform"],
                spectrum_id=x["title"],
                run=run,
                score=x["score"],
                precursor_mz=x["precursor_mz"],
                retention_time=x["rtinseconds"],
                source="PepNovo+",
                metadata={
                    "n_shift": x["n_shift"],
                    "c_shift": x["c_shift"],
                    "score_pepnovo": x["score_pepnovo"],
                    "scans": x["scans_y"]
                }
            ),
            axis=1
        ).tolist()
    
    psmlist = PSMList(psm_list=psm_list)
    return psmlist