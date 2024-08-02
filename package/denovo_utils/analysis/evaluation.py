import numpy as np
import pandas as pd
from .pandas_utils import (
    collapse_casanovo_score,
    get_spectralis_score,
    get_psm_type,
    amino_acid_converter,
    drop_charge
)
from ..parsers.converters import SpectralisParser
import os
from .mixture_models import assign_ggm_clusters
from tqdm import tqdm

def get_precision_coverage_df(df: pd.DataFrame, source: str, score_col: str, ground_truth_source="sage"):
    all_spectra = len(df[df.source==ground_truth_source])
    selection = df[df.source==source].sort_values(score_col)
    correctness_bool = selection["correct_prediction"].to_numpy()
    scores = selection[score_col].to_numpy()

    metrics = {
        "TP": [],
        "FP": [],
        "FN": [],
        "U": [],
        "U_mgf": []
    }
    
    for threshold in scores:
        bool_positive = scores > threshold
        bool_negative = scores < threshold

        
        tp = len(scores[np.logical_and(bool_positive, correctness_bool)])
        fp = len(scores[np.logical_and(bool_positive, ~correctness_bool)])
        fn = len(scores[np.logical_and(bool_negative, correctness_bool)])

        unpredicted = len(scores[bool_negative])
        unpredicted_mgf = all_spectra - len(scores[bool_negative])

        metrics["TP"].append(tp)
        metrics["FP"].append(fp)
        metrics["FN"].append(fn)
        metrics["U"].append(unpredicted)
        metrics["U_mgf"].append(unpredicted_mgf)

    metrics_df = pd.DataFrame(metrics)
    metrics_df["precision"] = metrics_df.apply(
        lambda x: x["TP"]/(x["TP"]+x["FP"]),
        axis=1
    )
    metrics_df["recall"] = metrics_df.apply(
        lambda x: x["TP"]/(x["TP"]+x["FN"]),
        axis=1
    )
    metrics_df["coverage_mgf"] = metrics_df.apply(
        lambda x: (x["TP"]+x["FP"])/(x["TP"]+x["FP"]+x["U_mgf"]),
        axis=1
    )
    metrics_df["coverage"] = metrics_df.apply(
        lambda x: (x["TP"]+x["FP"])/(x["TP"]+x["FP"]+x["U"]),
        axis=1
    )
    return metrics_df.dropna(0)


def calculate_accuracy(df, source, subset_spectra):
    df_subset = df[
        (df["source"]==source) &
        (df["spectrum_id"].isin(subset_spectra))
    ]
    if len(df_subset)==0:
        return None, 0
    precision = df_subset["correct_prediction"].sum() / len(df_subset)

    coverage = len(df_subset)/len(subset_spectra)

    return precision, coverage

def calculate_all_accuracy_metrics(df, engines, spectra_type_dict):
    metrics_accuracy = {}
    metrics_coverage = {}
    extra_predictions = {}


    for engine in engines:
        accuracy_accepted, coverage_accepted = calculate_accuracy(
            df=df,
            source=engine,
            subset_spectra=spectra_type_dict["target_accepted"]
        )
        accuracy_rejected, coverage_rejected = calculate_accuracy(
            df=df,
            source=engine,
            subset_spectra=spectra_type_dict["target_rejected"]
        )
        accuracy_decoy, coverage_decoy = calculate_accuracy(
            df=df,
            source=engine,
            subset_spectra=spectra_type_dict["decoy"]
        )
        metrics_accuracy[engine] = {
            "target_accepted": float(accuracy_accepted),
            "target_rejected": float(accuracy_rejected),
            "decoy": float(accuracy_decoy)
        }
        metrics_coverage[engine] = {
            "target_accepted": float(coverage_accepted),
            "target_rejected": float(coverage_rejected),
            "decoy": float(coverage_decoy)
        }

        subset = df[df.source==engine]
        extra_predictions_abs = subset["correct_prediction"].isna().sum() 
        extra_predictions_pct = extra_predictions_abs / len(subset)
        extra_predictions[engine] = {
            "pct": float(extra_predictions_pct),
            "absolute": int(extra_predictions_abs)
        }
    return metrics_accuracy, metrics_coverage, extra_predictions



def subset_spectra_on_psmtype(df, source_truth="sage"):
    df_db = df[df["source"]==source_truth]
    spectra_target_accepted = df_db.loc[df_db["psm_type"]=="target_accepted", "spectrum_id"].tolist()
    spectra_target_rejected = df_db.loc[df_db["psm_type"]=="target_rejected", "spectrum_id"].tolist()
    spectra_decoy = df_db.loc[df_db["psm_type"]=="decoy", "spectrum_id"].tolist()
    ground_truth = df_db.loc[:, ["peptide", "spectrum_id"]].set_index("spectrum_id").to_dict("index")

    return {
        "target_accepted": spectra_target_accepted,
        "target_rejected": spectra_target_rejected,
        "decoy": spectra_decoy,
        "ground_truth": ground_truth
    }

def evaluate_prediction(row, ground_truth):
    try:
        ground_truth_peptidoform = ground_truth[row["spectrum_id"]]["peptide"]
        return ground_truth_peptidoform == row["peptide"]

    except:
        return None

    
def load_and_preprocess(root, filename):
    mgf_path = os.path.join(root, "mgf_filtered", filename+".mgf")
    results_dir = os.path.join(root, "denovo_results")

    ### Parse raw output
    parser_spectralis = SpectralisParser(
        mgf_path=mgf_path,
        results_dir=results_dir
    )

    # Casanovo, instanovo, pepnet, contranovo ran together
    parser_spectralis.parse(
        path_spectralis=os.path.join(
            results_dir,
            "refinement/spectralis/pt1", filename + "_annotated_rescoring.csv"
        )
    )

    # NovoB, Novor, PepNovo+ ran together
    parser_spectralis.parse(
        path_spectralis=os.path.join(
            results_dir,
            "refinement/spectralis/pt2", filename + "_annotated_rescoring.csv"
        )
    )

    # Sage results ran separately
    parser_spectralis.parse(
        path_spectralis=os.path.join(
            results_dir,
            "refinement/spectralis/pt3", filename + "_annotated_rescoring.csv"
        )
    )

    ### Preprocess
    df = parser_spectralis.dataframe
    mapping_aa = {
        "I": "L",
        "N[UNIMOD:7]": "D",
        "Q[UNIMOD:7]": "E",
        "UNLMOD": "UNIMOD" # Otherwise the I in UNIMOD is replaced with L
    }

    df["spectralis_score"] = df.apply(get_spectralis_score, axis=1)
    df["score"] = df.apply(collapse_casanovo_score, axis=1)
    df["psm_type"] = df.apply(get_psm_type, axis=1)

    #df["peptidoform"] = df.apply(drop_charge, axis=1) wby would I do this?
    df["peptide"] = df.apply(
        lambda x: amino_acid_converter(x, mapping_aa),
        axis=1
    )
    df["peptide"] = df.apply(lambda x: x["peptide"].proforma, axis=1)

    # Prevent PSM-pydantic errors
    df["rank"] = 1
    df["run"] = filename

    return df

def annotate_GMM_clusters(df, db_engine="sage", only_df=True):
    
    new_df = []
    metrics_filter = {}
    df["gmm_cluster"] = None

    for engine in tqdm(df["source"].unique()):
        subset = df[df["source"]==engine].copy()
        

        if engine==db_engine:
            # Prevent filtering on the ground truth data
            # Should be filtered with standard TD x% FDR
            subset["gmm_cluster"] = 1
            new_df.append(subset)
            continue
        
        labels, mean, covs, weights = assign_ggm_clusters(
            subset["spectralis_score"].to_numpy(), n_clusters=2
        )

        subset["gmm_cluster"] = labels
        
        highest_cluster = subset.groupby("gmm_cluster")["spectralis_score"].max().idxmax()

        if highest_cluster==1:
            new_df.append(subset)
        else:
            subset["gmm_cluster"] = subset["gmm_cluster"].apply(
                lambda x: 1 if x==0 else 0
            )
            new_df.append(subset)

        if not only_df:
            cutoff = np.mean(
                [
                    subset.loc[subset["gmm_cluster"]==1, "spectralis_score"].min(),
                    subset.loc[subset["gmm_cluster"]==0, "spectralis_score"].max()
                ]
            )
            metrics_filter[engine] = {
                "spectra_all_n": len(subset),
                "spectra_accepted_n": len(subset[subset["gmm_cluster"]==1]),
                "spectra_filtered_n": len(subset[subset["gmm_cluster"]==0]),
                "spectralis_cutoff": float(cutoff),
                # "labels": [labels],
                "mean": [float(i) for i in mean],
                "covs": [float(i) for i in covs],
                "weights": [float(i) for i in weights]
            }
    
    new_df = pd.concat(new_df, ignore_index=True)
    if not only_df:
        return new_df, metrics_filter
    return new_df, None


def filter_gmm_cluster(df, keep_cluster_label=1):
    return df[df["gmm_cluster"]==keep_cluster_label].reset_index(drop=True)


def calculate_metrics(df, engines):
    db_result_dict = subset_spectra_on_psmtype(df)
    df["correct_prediction"] = df.apply(
        lambda x: evaluate_prediction(
            x,
            db_result_dict["ground_truth"]
        ),
        axis=1
    )
    metrics_accuracy, metrics_coverage, extra_predictions = calculate_all_accuracy_metrics(
        df,
        engines=engines,
        spectra_type_dict=db_result_dict
    )
    return {
        "accuracy": metrics_accuracy,
        "coverage": metrics_coverage,
        "extra_predictions": extra_predictions
    }

