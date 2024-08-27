"""Functions for calculating performance metrics of de novo search engines."""

import os
from typing import Optional

import numpy as np
import pandas as pd
from tqdm import tqdm

from ..parsers.converters import SpectralisParser
from ..utils.pandas import (
    amino_acid_converter,
    collapse_casanovo_score,
    evaluate_prediction,
    get_psm_type,
    get_spectralis_score,
)
from .mixture_models import assign_ggm_clusters


def get_precision_coverage_df(
    df: pd.DataFrame, source: str, score_col: str, 
    correctness_col: str,
    ground_truth_source: str = "sage"
) -> pd.DataFrame:
    """
    Calculate and return precision coverage dataframe for plotting.

    Parameters
    ----------
    df: pd.DataFrame
        A dataframe containing data for evaluation.
        Should contain the 'source' column with the values source and
        ground_truth_source in there, and the column 'score' with floats.
    source: str
        The search engine under evaluation.
        The value in column 'source' for which the metrics should be calculated.
    score_col: str
        The column name of the score column. Used to set score thresholds.
    ground_truth_source: str (default: sage)
        The name of the search engine (in source column) used as ground truth.

    Return
    ------
    pd.DataFrame
        The precision-coverage table.

    Notation
    --------
    - TP: True positives (PSMs above score theshold AND correct)
    - FP: False positives (PSMs above score threshold AND incorrect)
    - FN: False negatives (PSMs below score threshold AND correct)
    - U: All PSMs below the score threshold (PSM)
    - U_mgf: Unpredicted spectra in the complete set of spectra
    (error in calculation here??? Check!)

    Calculations
    ------------
    - precision: TP / (TP + FP)
    - recall: TP / (TP + FN)
    - coverage: (TP + FP) / (TP + FP + U)
    - coverage_mgf: (TP + FP) / (TP + FP + U_mgf)
    """
    all_spectra = len(df[df.source == ground_truth_source])
    selection = df[df.source == source].sort_values(score_col)
    correctness_bool = selection[correctness_col].to_numpy()
    scores = selection[score_col].to_numpy()

    metrics = {"TP": [], "FP": [], "FN": [], "U": [], "U_mgf": []}

    for threshold in scores:
        bool_positive = scores > threshold
        bool_negative = scores < threshold

        tp = len(scores[np.logical_and(bool_positive, correctness_bool)])
        fp = len(scores[np.logical_and(bool_positive, ~correctness_bool)])
        fn = len(scores[np.logical_and(bool_negative, correctness_bool)])

        unpredicted = len(scores[bool_negative])
        unpredicted_mgf = all_spectra - len(scores[bool_positive])

        metrics["TP"].append(tp)
        metrics["FP"].append(fp)
        metrics["FN"].append(fn)
        metrics["U"].append(unpredicted)
        metrics["U_mgf"].append(unpredicted_mgf)

    metrics_df = pd.DataFrame(metrics)
    metrics_df["precision"] = metrics_df.apply(
        lambda x: x["TP"] / (x["TP"] + x["FP"]), axis=1
    )
    metrics_df["recall"] = metrics_df.apply(
        lambda x: x["TP"] / (x["TP"] + x["FN"]), axis=1
    )
    metrics_df["coverage_mgf"] = metrics_df.apply(
        lambda x: (x["TP"] + x["FP"]) / (x["TP"] + x["FP"] + x["U_mgf"]), axis=1
    )
    metrics_df["coverage"] = metrics_df.apply(
        lambda x: (x["TP"] + x["FP"]) / (x["TP"] + x["FP"] + x["U"]), axis=1
    )
    return metrics_df

def get_threshold(df, threshold=.9):

    for i, precision in enumerate(df.precision.to_numpy()):
        if precision > threshold:
            return i
        
def profile_accuracy_thresholds(
    df,
    engine,
    score_col="score",
    thresholds = [.99, .95, .9, .85, .8]
):
    df_subset = df[
        df.source==engine
    ].sort_values(score_col)
    
    df["correct_prediction"] = df["match_type"].apply(lambda x: x=="Match")
    df["correct_prediction_lenient"] = df["match_type"].apply(lambda x: x!="Worse")
    metrics_score = get_precision_coverage_df(
        df,
        source=engine,
        score_col=score_col,
        correctness_col="correct_prediction",
        ground_truth_source="sage"
    )
    metrics_lenient = get_precision_coverage_df(
        df,
        source=engine,
        score_col=score_col,
        correctness_col="correct_prediction_lenient",
        ground_truth_source="sage"
    )
    
    score_profile = {"Accuracy": {}, "Accuracy_lenient": {}}
    for threshold in thresholds:
        threshold_profile = {}
        threshold_profile_lenient = {}

        iloc_threshold = get_threshold(df=metrics_score, threshold=threshold)
        iloc_threshold_lenient = get_threshold(df=metrics_lenient, threshold=threshold)

        threshold_profile["recall"] = metrics_score.iloc[iloc_threshold]["recall"]
        threshold_profile["coverage"] = metrics_score.iloc[iloc_threshold]["coverage"]
        threshold_profile["score"] = df_subset.iloc[iloc_threshold][score_col]
        
        threshold_profile_lenient["recall"] = metrics_lenient.iloc[iloc_threshold_lenient]["recall"]
        threshold_profile_lenient["coverage"] = metrics_lenient.iloc[iloc_threshold_lenient]["coverage"]
        threshold_profile_lenient["score"] = df_subset.iloc[iloc_threshold_lenient][score_col]
        partition = df_subset.iloc[iloc_threshold_lenient:]["match_type"].value_counts(normalize=True)
        for match_type in ["Match", "Better", "Worse", "Isobaric"]:
            try:
                threshold_profile_lenient[match_type] = partition[match_type]
            except KeyError:
                threshold_profile_lenient[match_type] = 0.0
    
        score_profile["Accuracy"][threshold] = threshold_profile
        score_profile["Accuracy_lenient"][threshold] = threshold_profile_lenient
    return score_profile

def calculate_accuracy(
    df: pd.DataFrame, source: str, subset_spectra: list
) -> tuple[float, float]:
    """
    Calculate the accuracy on a subset of spectra.

    Requires columns 'source', 'spectrum_id', and 'correct_prediction'.
    Filters the dataframe df on value set as source and the list of spectra
    in subset_spectra.

    Parameters
    ----------
    df: pd.DataFrame
        Dataframe containing columns as specified above.
    source: str
        The name of the search engine for which to compute accuracy.
    subset_spectra: list
        A list of spectra that should be included for accuracy calculation.
        This is most probably extracted based on the ground truth set of spectra.

    Return
    ------
    precision: float
        The precision calculated as: #correct_predictions / #predictions
    coverage: float
        The coverage calculated as #predictions / #spectra
    """
    df_subset = df[(df["source"] == source) & (df["spectrum_id"].isin(subset_spectra))]
    if len(df_subset) == 0:
        return None, 0
    precision = df_subset["correct_prediction"].sum() / len(df_subset)

    coverage = len(df_subset) / len(subset_spectra)

    return precision, coverage


def calculate_all_accuracy_metrics(
    df: pd.DataFrame, engines: list, spectra_type_dict: dict
) -> tuple[dict, dict, dict]:
    """
    Calculate accuracy metrics for subsets target_accepted, rejected and decoys.

    df requires columns 'source', 'spectrum_id', and 'correct_prediction'.
    Filters the dataframe df on value set as source and the list of spectra
    in subset_spectra.

    Parameters
    ----------
    df: pd.DataFrame
        Dataframe containing columns as specified above
    engines: list
        List of engines (values in source column).
    spectra_type_dict: dict
        Dictionary containing keys 'target_rejected', 'target_accepted', and 'decoy'.
        Values should be a list of spectrum_ids that were matched as specified by keys.

    Return
    ------
    metrics_accuracy: dict
        Accuracy score for each PSM-type.
    metrics_coverage: dict
        Coverage score for each PSM-type.
    extra_predictions: dict
        'pct' and 'absolute' relate to number of spectra predicted not matched with a
        peptide sequence with the database search (ground truth).
    """
    metrics_accuracy = {}
    metrics_coverage = {}
    extra_predictions = {}

    for engine in engines:
        accuracy_accepted, coverage_accepted = calculate_accuracy(
            df=df, source=engine, subset_spectra=spectra_type_dict["target_accepted"]
        )
        accuracy_rejected, coverage_rejected = calculate_accuracy(
            df=df, source=engine, subset_spectra=spectra_type_dict["target_rejected"]
        )
        accuracy_decoy, coverage_decoy = calculate_accuracy(
            df=df, source=engine, subset_spectra=spectra_type_dict["decoy"]
        )
        metrics_accuracy[engine] = {
            "target_accepted": float(accuracy_accepted),
            "target_rejected": float(accuracy_rejected),
            "decoy": float(accuracy_decoy),
        }
        metrics_coverage[engine] = {
            "target_accepted": float(coverage_accepted),
            "target_rejected": float(coverage_rejected),
            "decoy": float(coverage_decoy),
        }

        subset = df[df.source == engine]
        extra_predictions_abs = subset["correct_prediction"].isna().sum()
        extra_predictions_pct = extra_predictions_abs / len(subset)
        extra_predictions[engine] = {
            "pct": float(extra_predictions_pct),
            "absolute": int(extra_predictions_abs),
        }
    return metrics_accuracy, metrics_coverage, extra_predictions


def subset_spectra_on_psmtype(df: pd.DataFrame, source_truth: str = "sage") -> dict:
    """
    Create dictionary for each psm-type with spectrum ids.

    Parameters
    ----------
    df: pd.DataFrame
        DataFrame with columns 'source', 'psm_type', 'peptide', and 'spectrum_id'
    source_truth: str
        The value in column 'source' corresponding with ground truth.

    Return
    ------
    dict:
        A dictionary with a list of spectrum_ids for each psm_type.
        Only the ground_truth value is a dictionary.
    """
    df_db = df[df["source"] == source_truth]
    spectra_target_accepted = df_db.loc[
        df_db["psm_type"] == "target_accepted", "spectrum_id"
    ].tolist()
    spectra_target_rejected = df_db.loc[
        df_db["psm_type"] == "target_rejected", "spectrum_id"
    ].tolist()
    spectra_decoy = df_db.loc[df_db["psm_type"] == "decoy", "spectrum_id"].tolist()
    ground_truth = (
        df_db.loc[:, ["peptide", "spectrum_id"]]
        .set_index("spectrum_id")
        .to_dict("index")
    )

    return {
        "target_accepted": spectra_target_accepted,
        "target_rejected": spectra_target_rejected,
        "decoy": spectra_decoy,
        "ground_truth": ground_truth,
    }


def load_and_preprocess(root: str, filename: str) -> pd.DataFrame:
    """Crap function with hard-coded stuff..."""
    mgf_path = os.path.join(root, "mgf_filtered", filename + ".mgf")
    results_dir = os.path.join(root, "denovo_results")

    # 1. Parse raw output
    parser_spectralis = SpectralisParser(mgf_path=mgf_path, results_dir=results_dir)

    # Casanovo, instanovo, pepnet, contranovo ran together
    parser_spectralis.parse(
        path_spectralis=os.path.join(
            results_dir,
            "refinement/spectralis/pt1",
            filename + "_annotated_rescoring.csv",
        )
    )

    # NovoB, Novor, PepNovo+ ran together
    parser_spectralis.parse(
        path_spectralis=os.path.join(
            results_dir,
            "refinement/spectralis/pt2",
            filename + "_annotated_rescoring.csv",
        )
    )

    # Sage results ran separately
    parser_spectralis.parse(
        path_spectralis=os.path.join(
            results_dir,
            "refinement/spectralis/pt3",
            filename + "_annotated_rescoring.csv",
        )
    )

    # 2. Preprocess
    df = parser_spectralis.dataframe
    mapping_aa = {
        "N[UNIMOD:7]": "D",
        "Q[UNIMOD:7]": "E",
        "I": "L",
        "UNLMOD": "UNIMOD",  # Otherwise the I in UNIMOD is replaced with L
    }

    df["spectralis_score"] = df.apply(get_spectralis_score, axis=1)
    df["score"] = df.apply(collapse_casanovo_score, axis=1)
    df["psm_type"] = df.apply(get_psm_type, axis=1)

    # df["peptidoform"] = df.apply(drop_charge, axis=1) wby would I do this?
    df["peptide"] = df.apply(lambda x: amino_acid_converter(x, mapping_aa), axis=1)
    df["peptide"] = df.apply(lambda x: x["peptide"].proforma, axis=1)

    # Prevent PSM-pydantic errors
    df["rank"] = 1
    df["run"] = filename

    return df


def annotate_GMM_clusters(
    df: pd.DataFrame, score_col: str, db_engine: str = "sage", only_df: bool = True
) -> tuple[pd.DataFrame, Optional[dict]]:
    """
    REWRITE THIS FUNCTION.

    Annotate a set of PSMs dependent on allocation to a model defined
    with Gaussian Mixture Models.

    Parameters
    ----------
    df: pd.DataFame
        Dataframe containing columns 'source', score_col.
        Source defines the engine of which its predictions must be modelled,
        based on the scores in score_col.
    score_col: str
        The column name of the score column.
    db_engine: str (default: sage)
        The value in the 'source' column defining the engine of which the
        PSMs should be modelled with GMM.
    only_df: bool (default: True)
        Whether to only return the dataframe with GMM labels, or to also
        include a dictionary with a bunch of GMM-specific metrics.

    Returns
    -------
    df: pd.DataFrame
        DataFrame with extra gmm_cluster where 1 is an accepted and 0 a rejected PSM.
    metrics: dict
        Keys include:
        - spectra_all_n
        - spectra_accepted_n
        - spectra_filtered_n
        - score_cutoff
        - mean
        - covs
        - weights
    """
    new_df = []
    metrics_filter = {}
    df["gmm_cluster"] = None

    for engine in tqdm(df["source"].unique()):
        subset = df[df["source"] == engine].copy()

        if engine == db_engine:
            # Prevent filtering on the ground truth data
            # Should be filtered with standard TD x% FDR
            subset["gmm_cluster"] = 1
            new_df.append(subset)
            continue

        labels, mean, covs, weights = assign_ggm_clusters(
            subset[score_col].to_numpy(), n_clusters=2
        )

        subset["gmm_cluster"] = labels

        highest_cluster = subset.groupby("gmm_cluster")[score_col].max().idxmax()

        if highest_cluster == 1:
            new_df.append(subset)
        else:
            subset["gmm_cluster"] = subset["gmm_cluster"].apply(
                lambda x: 1 if x == 0 else 0
            )
            new_df.append(subset)

        if not only_df:
            cutoff = np.mean(
                [
                    subset.loc[subset["gmm_cluster"] == 1, score_col].min(),
                    subset.loc[subset["gmm_cluster"] == 0, score_col].max(),
                ]
            )
            metrics_filter[engine] = {
                "spectra_all_n": len(subset),
                "spectra_accepted_n": len(subset[subset["gmm_cluster"] == 1]),
                "spectra_filtered_n": len(subset[subset["gmm_cluster"] == 0]),
                "score_cutoff": float(cutoff),
                # "labels": [labels],
                "mean": [float(i) for i in mean],
                "covs": [float(i) for i in covs],
                "weights": [float(i) for i in weights],
            }

    new_df = pd.concat(new_df, ignore_index=True)
    if not only_df:
        return new_df, metrics_filter
    return new_df, None


def filter_gmm_cluster(df: pd.DataFrame, keep_cluster_label: int = 1) -> pd.DataFrame:
    """Filter and return the dataframe by gmm_cluster==keep_cluster_label."""
    return df[df["gmm_cluster"] == keep_cluster_label].reset_index(drop=True)


def calculate_metrics(df: pd.DataFrame, engines: list[str]) -> dict:
    """Calculate accuracy metrics for the specified set of engines."""
    db_result_dict = subset_spectra_on_psmtype(df)
    df["correct_prediction"] = df.apply(
        lambda x: evaluate_prediction(x, db_result_dict["ground_truth"]), axis=1
    )
    metrics_accuracy, metrics_coverage, extra_predictions = (
        calculate_all_accuracy_metrics(
            df, engines=engines, spectra_type_dict=db_result_dict
        )
    )
    return {
        "accuracy": metrics_accuracy,
        "coverage": metrics_coverage,
        "extra_predictions": extra_predictions,
    }
