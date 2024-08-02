import pandas as pd
from psm_utils import Peptidoform
from psm_utils.peptidoform import PeptidoformException
import logging
import re
import pyopenms as oms

def parse_peptidoform(peptide: str, mapping: dict, max_length=30):
    peptide_parsed = peptide
    for k, v in mapping.items():
        if ("-" in v) and (not peptide_parsed.startswith(k)):
            peptide_parsed = peptide_parsed.replace(k, v[:-1])
        else:
            peptide_parsed = peptide_parsed.replace(k, v)

    try:
        peptidoform = Peptidoform(peptide_parsed)
        if (len(peptidoform) > max_length) or (peptidoform.precursor_charge > 6) or (len(peptidoform) < 2):
            return None
        return peptidoform
    except:
        logging.warning(f"Failed to parse: {peptide}")
        return None
