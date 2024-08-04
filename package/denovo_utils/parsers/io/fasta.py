import re

import pandas as pd
from Bio import SeqIO

from ...utils.pandas import row_to_seqrecord


class FastaHandler:
    def __init__(self):
        self.read_file = False

    def read(self, path_fasta) -> None:
        entries = []
        seq_record = SeqIO.parse(path_fasta, format="fasta")
        for entry in seq_record:
            entry_dict = self._parse_fasta_header(entry.description)
            entry_dict["sequence"] = entry.seq

            if not entry_dict["protein_id_full"] == entry.id:
                raise Exception(entry_dict)
            entries.append(entry_dict)
        self.dataframe = pd.DataFrame(entries)
        self.read_file = True

    def to_dict(self, key, value):
        """
        Convert a fasta file to a dictionary based on a key and value column.

        Parameters
        ----------
        key: str
            The dataframe column name used as the dictionary key
        value: str
            The dataframe column name used as the dictionary values
        """
        if not self.read_file:
            raise Exception("Read a fasta file first with the 'read' method")
        if key not in self.dataframe.columns or value not in self.dataframe.columns:
            raise Exception(
                "Make sure the key and value parameters are a valid column name. Columns: {}".format(
                    self.dataframe.columns
                )
            )

        return self.dataframe[[key, value]].groupby(key)[value].apply(list).to_dict()

    def _parse_fasta_header(self, fasta_header):
        """
        Parse a fasta header into a dictionary.

        Parsed elements include:
            'name'
            'protein_id'
            'protein_id_full
            'protein_description'
            'organism'
            'gene_name'
            'protein_existence'
            'sequence_version'
        """
        # Regex patterns for parsing
        entry_name_pattern = r"^(tr|sp|\w+)\|(\w+)\|"
        protein_id_full_pattern = r"^(tr|sp|\w+\|[\w_]+(\|\S*)?)(?=\s|$)"
        description_pattern = r" (.*?) (\w+=)"
        tag_pattern = r"(\w+)=([^=]+?)(?=\s\w+=|$)"

        # Dictionary for mapped entries
        entry_name_mapping = {"tr": "TrEMBL", "sp": "SwissProt"}

        # Parse entry name and protein ID
        entry_name_match = re.search(entry_name_pattern, fasta_header)
        entry_name = entry_name_match.group(1) if entry_name_match else "Other"
        entry_name = entry_name_mapping.get(entry_name, "Other")
        protein_id = entry_name_match.group(2) if entry_name_match else None

        # Parse full protein ID
        protein_id_full_match = re.search(protein_id_full_pattern, fasta_header)
        protein_id_full = (
            protein_id_full_match.group(1) if protein_id_full_match else None
        )

        # Parse protein description
        description_match = re.search(description_pattern, fasta_header)
        protein_description = description_match.group(1) if description_match else None

        # Parse optional elements
        tags = re.findall(tag_pattern, fasta_header)
        tag_dict = {tag: value for tag, value in tags}

        organism = tag_dict.get("OS")
        gene_name = tag_dict.get("GN")
        protein_existence = tag_dict.get("PE")
        sequence_version = tag_dict.get("SV")

        # If the regex of the protein name is not in UniProt format
        if isinstance(protein_id_full, type(None)):
            if len(fasta_header.split()) == 1:
                protein_id_full = fasta_header
            elif isinstance(protein_description, type(None)):
                fasta_list = fasta_header.split()
                protein_id_full = fasta_list[0]
                protein_description = " ".join(fasta_list[1:])
            else:
                raise Exception(f"{fasta_header} could not be parsed.")

        # Construct result dictionary
        result = {
            "name": entry_name,
            "protein_id": protein_id,
            "protein_id_full": protein_id_full,
            "protein_description": protein_description,
            "organism": organism,
            "gene_name": gene_name,
            "protein_existence": protein_existence,
            "sequence_version": sequence_version,
        }
        return result

    @property
    def organisms(self):
        organisms = self.dataframe["organism"].to_list()
        if self.dataframe["organism"].isna().sum() > 0:
            print(
                "Warning: Some sequences have no annotation for organism. Please fix."
            )
        return organisms

    def write(self, boolean_filter: list[bool], out_path: str):
        df_out = self.dataframe.loc[boolean_filter]
        seq_records = df_out.apply(row_to_seqrecord, axis=1)

        with open(out_path, "w") as f:
            _ = SeqIO.write(sequences=seq_records, handle=f, format="fasta")


class FastaWriter:
    def __init__(self):
        pass
