from psm_utils.io.mzid import MzidReader
from psm_utils.io import read_file
from psm_utils import Peptidoform
from pathlib import Path
from typing import Union
import os

class MSGFReader(MzidReader):
    def __init__(self, filename: str | Path, *args, score_key: str = None, **kwargs) -> None:
        super().__init__(filename, *args, **kwargs)

    @staticmethod
    def _parse_peptidoform(seq: str, modification_list: list[dict], charge: Union[int, None]):
        """Parse mzid sequence and modifications to Peptidoform."""
        peptide = [""] + list(seq) + [""]

        # Add modification labels
        for mod in modification_list:
            try:
                peptide[int(mod["location"])] += f"[{mod['name']}]"
            except:
                if mod['unknown modification'] == 'combination':

                    peptide[int(mod["location"])] += f"[+{mod['monoisotopicMassDelta']}]"
                else:
                    raise Exception(f'Other bug related to {modification_list}')

        # Add dashes between residues and termini, and join sequence
        peptide[0] = peptide[0] + "-" if peptide[0] else ""
        peptide[-1] = "-" + peptide[-1] if peptide[-1] else ""
        proforma_seq = "".join(peptide)

        # Add charge state
        if charge:
            proforma_seq += f"/{charge}"

        return Peptidoform(proforma_seq)
    
def msgf_parser_mzid(
        result_path: str,
):
    result_path = os.path.splitext(result_path)[0] + '.mzid'
    psmlist = MSGFReader(
        result_path
    ).read_file()

    # Fix retention time
    psmlist['retention_time'] = psmlist['retention_time'] / 60
    return psmlist

def msgf_parser_parquet(
        result_path: str,
):
    result_path = os.path.splitext(result_path)[0] + '.parquet'
    
    psmlist = read_file(
        result_path,
        filetype='parquet'
    )
    # Fix retention time
    psmlist['retention_time'] = psmlist['retention_time'] / 60
    return psmlist

def msgf_parser(
        result_path: str,
        mgf_path: str,
        mapping: dict,
        max_length=30,
        **kwargs
):
    if result_path.endswith('.mzid'):
        return msgf_parser_mzid(
            result_path=result_path
        )

    elif result_path.endswith('.parquet'):
        return msgf_parser_parquet(
            result_path=result_path
        )

    else:
        raise Exception(f"Unknown extension for MSGF parsing: {result_path}")