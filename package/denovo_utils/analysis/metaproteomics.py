import pandas as pd
import seaborn as sns
from ..parsers import DenovoEngineConverter
from .pandas_utils import get_psm_type

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

def count_analysis(path_file: str, psm_types: list, fasta_dict: dict, engine="sage", plot=True):
    """
    Perform a species-specific peptide analysis.

    Goal of the analysis is to check score distributions of PSMs of different species.

    Parameters
    ----------
    path_file: str
        The path to the sage search results
    psm_types: list[str]
        list with 'decoy', 'target_rejected', 'target_accepted', dependent on which subset
        should be included
    fasta_dict: dict
        A organism-protein mapping with keys (organism) and values (protein ids).
        Can be generated with the FastaReader class.

    Returns
    -------
    pd.DataFrame
        Dataframe containing the search results, filtered as specified
    pd.DataFrame
        A count table of species-specific peptides
    """
    # Read sage results
    parser = DenovoEngineConverter.select(label=engine)
    df_sage = parser.parse(
        result_path=path_file,
        mgf_path=""
    ).to_dataframe()

    # Filter on target_accepted
    df_sage["psm_type"] = df_sage.apply(get_psm_type, axis=1)
    df_sage_filtered = df_sage[df_sage["psm_type"].isin(psm_types)]

    # Get species specific peptides
    df_sage_filtered["species"] = df_sage_filtered.progress_apply(
        lambda x: get_species_specificity(x["protein_list"], fasta_dict),
        axis=1
    )
    df_sage_filtered["species_n"] = df_sage_filtered["species"].apply(len)

    # Count species specificity
    count_species_specific = count_species_peptides(df_sage_filtered)
    count_species = count_species_peptides(
        df=df_sage_filtered,
        species_list=count_species_specific["species"].tolist(),
        specific=False
    )
    count_table = pd.merge(count_species_specific, count_species)

    # Plot score distribution by species
    df_sage_unique_species = df_sage_filtered[df_sage_filtered["species_n"]==1].copy()
    df_sage_unique_species["species"] = df_sage_unique_species["species"].apply(lambda x: x[0])

    if plot:
        sns_plot = sns.kdeplot(
            df_sage_unique_species,
            x="score",
            hue="species",
            common_norm=False
        )
        sns.move_legend(sns_plot, "upper left", bbox_to_anchor=(1, 1))
    return df_sage_filtered, count_table