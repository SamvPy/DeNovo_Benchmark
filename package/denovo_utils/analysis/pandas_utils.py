from psm_utils import Peptidoform

def get_decoy_status(row):
    proteins = row["protein_list"]
    for protein in proteins:
        if "decoy" in protein:
            return True
    return False

def collapse_casanovo_score(row):
    if (row["source"] in ["Casanovo4.2.0", "ContraNovo"]) & (row["score"]<0):
        return 1+row["score"]
    return row["score"]

def get_spectralis_score(row):
    try:
        return row["rescoring_features"]["spectralis_score"]
    except:
        return None
    
def get_psm_type(row, cutoff=.01):
    if not isinstance(row["is_decoy"], bool):
        return None
    elif row["is_decoy"]:
        return "decoy"
    elif row["qvalue"]<cutoff:
        return "target_accepted"
    return "target_rejected"

def amino_acid_converter(row, mapping):
    peptidoform = row["peptidoform"].proforma
    for k, v in mapping.items():
        peptidoform = peptidoform.replace(k,v)
    return Peptidoform(peptidoform)

def drop_charge(row):
    return Peptidoform(row["peptidoform"].proforma.split("/")[0])