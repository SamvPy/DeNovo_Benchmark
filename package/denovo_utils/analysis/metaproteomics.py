import pandas as pd
import seaborn as sns
from ..parsers import DenovoEngineConverter
from ..utils.pandas import get_psm_type, get_species_specificity, count_species_peptides

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