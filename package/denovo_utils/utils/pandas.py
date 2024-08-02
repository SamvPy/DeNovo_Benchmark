from psm_utils import Peptidoform
import pandas as pd

def get_decoy_status(row, decoy_strings=["DECOY", "rev"]):
    proteins = row["protein_list"]
    for protein in proteins:
        for decoy_string in decoy_strings:
            if decoy_string in protein:
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

def get_species_specificity(protein_list, fasta_dict):
    species = []
    for protein in protein_list:
        for organism, org_spec_proteins in fasta_dict.items():
            if organism in species:
                continue
            if protein in org_spec_proteins:
                species.append(organism)
    return species

def count_species_peptides(df: pd.DataFrame, species_list=[], specific=True) -> pd.DataFrame:
    if specific:
        count_table = df.loc[df["species_n"]==1, "species"].value_counts(0).reset_index().rename(
            columns={"species": "count_specific", "index": "species"}
        )
        count_table["species"] = count_table["species"].apply(lambda x: x[0])
        return count_table
    
    count_dict = {
        species: 0 for species in species_list
    }
    def count_species(species_list, count_dict):
        for species in species_list:
            count_dict[species] += 1
    
    df["species"].progress_apply(lambda x: count_species(x, count_dict=count_dict))
    return pd.DataFrame({k: [v] for k, v in count_dict.items()}).melt().rename(columns={
        "variable": "species", "value": "count_all"
    })