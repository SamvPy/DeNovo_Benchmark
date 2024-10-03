"""Utilities useable on rows in a pandas dataframe with apply method."""

import re
from typing import Literal, Optional

import pandas as pd
from Bio.SeqRecord import SeqRecord
from psm_utils import Peptidoform

# SCORING-RELATED UTILITIES ###


def collapse_casanovo_score(row: pd.Series) -> int:
    """
    Ignore precursor mass shift tag (-1) of the casanovo score.

    Essentially, moves negative scores to the positive range.
    This happens when the engine tags PSMs as those
    not fitting with the precursor mass.

    Uses 'source' column to check for the search engine to parse scores for.
    Only parses 'source' with values 'Casanovo4.2.0' or 'ContraNovo'.

    Uses 'score' column as column to change if the value is negative.

    Parameter
    ---------
    row: pd.Series
        A row in a dataframe, or a series.

    Return
    ------
    int:
        1 + score if negative, otherwise original score.
    """
    if (row["source"] in ["Casanovo4.2.0", "ContraNovo"]) & (row["score"] < 0):
        return 1 + row["score"]
    return row["score"]


def get_spectralis_score(row: pd.Series) -> Optional[int]:
    """
    Return the spectralis score from the rescoring_features dict in a psm_list.

    Parameter
    ---------
    row: pd.Series
        A row in a dataframe, or a series.

    Return
    ------
    int, None
        The spectralis score if in the rescoring_features dictionary, otherwise None.
    """
    try:
        return row["rescoring_features"]["spectralis_score"]
    except KeyError:
        return None


# SEARCH RESULT RELATED UTILITIES ###


def get_decoy_status(
    row: pd.Series, decoy_strings: list[str] = ["DECOY", "rev"]
) -> bool:
    """
    Return decoy status of a PSM as boolean based on the protein identifier.

    Uses 'protein_list' column, containing a list of protein ids.

    Parameter
    ---------
    row: pd.Series
        A row in a dataframe, or a series.
    decoy_strings: list[str] (default: ['DECOY', 'rev'])
        String subsets present in protein identifiers indicating decoy status.

    Return
    ------
    bool
        True if decoy, otherwise False
    """
    if not isinstance(row["protein_list"], list):
        return None
    proteins = row["protein_list"]
    for protein in proteins:
        for decoy_string in decoy_strings:
            if decoy_string in protein:
                return True
    return False


def get_psm_type(
    row: pd.Series, cutoff: int = 0.01
) -> Literal["decoy", "target_rejected", "target_accepted"]:
    """
    Return the PSM-type based on 'is_decoy' and 'qvalue' columns.

    Parameter
    ---------
    row: pd.Series
        A row in a dataframe, or a series.
    cutoff: int (default: 0.01)
        Q-value cutoff value indicating the FDR-threshold.

    Return
    ------
    str
        'decoy', 'target_accepted', or 'target_rejected'.
    """
    if not isinstance(row["is_decoy"], bool):
        return None
    elif row["is_decoy"]:
        return "decoy"
    elif row["qvalue"] < cutoff:
        return "target_accepted"
    return "target_rejected"


# OTHER UTILTIES


def amino_acid_converter(row: pd.Series, mapping: dict) -> Peptidoform:
    """
    Parse amino acids or substrings from a proforma peptide string.

    Uses the 'peptidoform' column.

    Parameter
    ---------
    row: pd.Series
        A row in a dataframe, or a series.
    mapping: dict
        Mapping of modifications or residues.

    Return
    ------
    psm_utils.Peptidoform
        Parsed peptide sequence in Peptidoform format.
    """
    peptidoform = row["peptidoform"].proforma
    for k, v in mapping.items():
        peptidoform = peptidoform.replace(k, v)
    return Peptidoform(peptidoform)


def drop_charge(row: pd.Series) -> Peptidoform:
    """
    Drop the charge from a peptidoform peptide sequence.

    Uses the 'peptidoform' column.

    Parameter
    ---------
    row: pd.Series
        A row in a dataframe, or a series.
    """
    return Peptidoform(row["peptidoform"].proforma.split("/")[0])


# FASTA-RELATED UTILITIES ###


def get_species_specificity(protein_list: list[str], fasta_dict: dict) -> list[str]:
    """
    Return the species matching with the identified peptide.

    Parameter
    ---------
    protein_list: list
        List of protein ids, potentially matching with the value list in
        fasta_dict.
    fasta_dict: dict
        A dictionary of species - protein_id mappings.
        Generatable with FastaHandler.to_dict().

    Return
    ------
    list[str]:
        The species that match with the identified peptide.
    """
    species = []
    for protein in protein_list:
        for organism, org_spec_proteins in fasta_dict.items():
            if organism in species:
                continue
            if protein in org_spec_proteins:
                species.append(organism)
    return species


def count_species_peptides(
    df: pd.DataFrame, species_list=[], specific=True
) -> pd.DataFrame:
    """
    Count the peptides matching with a species.

    Uses columns:
    - species_n (if specific=True): The number of species matching with a peptide id.
    - species: The list of matching species

    Parameters
    ----------
    df: pd.DataFrame
        The dataframe used to count species specific peptides.
    species_list: list[str]
        The species for which peptide-specific counts needs to be computed.
    specific: bool
        Whether to filter on proteotypicity between species.
    """
    if specific:
        count_table = (
            df.loc[df["species_n"] == 1, "species"]
            .value_counts()
            .reset_index()
            .rename(columns={"count": "count_specific", "index": "species"})
        )
        count_table["species"] = count_table["species"].apply(lambda x: x[0])
        return count_table

    count_dict = {species: 0 for species in species_list}

    def count_species(species_list, count_dict):
        for species in species_list:
            count_dict[species] += 1

    df["species"].progress_apply(lambda x: count_species(x, count_dict=count_dict))
    return (
        pd.DataFrame({k: [v] for k, v in count_dict.items()})
        .melt()
        .rename(columns={"variable": "species", "value": "count_all"})
    )


def row_to_seqrecord(row: pd.Series) -> SeqRecord:
    """
    Convert a parsed fasta entry to a biopython SeqRecord object.

    Makes it trivial to write fasta files with a list of SeqRecords.

    Parameter
    ---------
    row: pd.Series
        A row in a dataframe, or a series.

    Return
    ------
    Bio.SeqRecord.SeqRecord
        A SeqRecord object of the fasta entry.
    """
    prefix_mapping = {
        "protein_id_full": "",
        "protein_description": "",
        "organism": "OS=",
        "gene_name": "GN=",
        "protein_existence": "PE=",
        "sequence_version": "SV=",
    }

    sequence = row["sequence"]
    id_ = row["protein_id_full"]
    name = row["protein_id_full"]
    description_list = []
    for col in [
        "protein_id_full",
        "protein_description",
        "organism",
        "gene_name",
        "protein_existence",
        "sequence_version",
    ]:
        value = row[col]
        if value:
            description_list.append(prefix_mapping[col] + row[col])
    description = " ".join(description_list)
    return SeqRecord(seq=sequence, id=id_, name=name, description=description)


def prediction_to_seqrecord(row: pd.Series) -> SeqRecord:
    """
    Parse a row representing a de novo prediction to a biopython SeqRecord.

    The row should have columns 'peptidoform', 'source', and 'spectrum_id'.
    The peptidoform can have type str (proforma format) or Peptidoform.

    The OS= tag is annotated with whatever is in 'source'. The name of the entry
    is denovo_<scan_number>. Lastly, the protein is constructed as 
    dn|<sequence>|.

    Parameter
    ---------
    row: pd.Series
        A row in a dataframe, or a series.
    
    Return
    ------
    Bio.SeqRecord.SeqRecord
        The de novo prediction parsed as a SeqRecord object.
    """
    pattern_scan = r"scan=(\d+)"

    sequence = Peptidoform(row["peptidoform"]).sequence
    id_ = f"dn|{sequence}|"
    organism = "OS=" + row["source"]

    # Either just a number (indicating the scan number) or the complete spectrum_id
    try:
        scan_number = re.search(pattern_scan, row["spectrum_id"]).group(1)
    except:
        scan_number = row["spectrum_id"]
    name = f"denovo_{scan_number}"

    seqrecord = SeqRecord(
        seq=sequence,
        id=id_,
        name=name,
        description=" ".join([id_, name, organism])
    )
    return seqrecord



# Evaluation-related utilities


def evaluate_prediction(row: pd.Series, ground_truth: pd.DataFrame) -> Optional[bool]:
    """
    Evaluate the truthness of a de novo prediction by comparing it with ground truth.

    Parameter
    ---------
    row: pd.Series
        The row containing the spectrum_id and peptide prediction in
        columns 'spectrum_id' and 'peptide' respectively.
    ground_truth: pd.DataFrame
        A dataframe with columns 'spectrum_id' and 'peptide' which are
        used to extract the spectrum_id of the PSM and the peptide sequence.

    Return
    ------
    Optional[bool]
        True or False dependent on correctness, None if the spectrum was not given
        a PSM with the ground truth search.
    """
    try:
        ground_truth_peptidoform = ground_truth.loc[row["spectrum_id"], "peptide"]
        return ground_truth_peptidoform == row["peptide"]

    except Exception:
        return None
