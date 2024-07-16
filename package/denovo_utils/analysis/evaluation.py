import numpy as np
import pandas as pd

def get_precision_metrics(df: pd.DataFrame, source: str, score_col: str, ground_truth_source="sage"):
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