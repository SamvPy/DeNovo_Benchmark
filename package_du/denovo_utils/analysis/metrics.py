"""Methods to evaluate peptide-spectrum predictions. Borrowed from
Casanovo https://github.com/Noble-Lab/casanovo"""

"""Taken from denovo_benchmarks from PominovaMS (https://github.com/PominovaMS/denovo_benchmarks)"""

import re
import numpy as np
import pandas as pd
from typing import Dict, Iterable, List, Tuple
from ..parsers.constants import AA_MASSES
from sklearn.metrics import precision_recall_curve
import seaborn as sns
import matplotlib.pyplot as plt

# Method is borrowed from 'spectrum_utils' package
# (temporary removed 'spectrum_utils' dependence
# due to numba cache issues when running on HPC cluster).
def mass_diff(mz1, mz2, mode_is_da):
    """
    Calculate the mass difference(s).

    Parameters
    ----------
    mz1
        First m/z value(s).
    mz2
        Second m/z value(s).
    mode_is_da : bool
        Mass difference in Dalton (True) or in ppm (False).

    Returns
    -------
        The mass difference(s) between the given m/z values.
    """
    return mz1 - mz2 if mode_is_da else (mz1 - mz2) / mz2 * 10**6


def get_token_mass(
    token: tuple
) -> float:
    """
    Convert the amino acid to a mass while considering modifications as well.
    """
    aa, mods = token
    mass = AA_MASSES[aa]
    for mod in mods:
        if mod is None:
            continue
        mass += mod.mass
    return mass


def aa_match_prefix(
    peptide1: List[str],
    peptide2: List[str],
    aa_dict: Dict[str, float] = {},
    cum_mass_threshold: float = 50,
    ind_mass_threshold: float = 20,
) -> Tuple[np.ndarray, bool, Tuple[np.ndarray]]:
    """
    Find the matching prefix amino acids between two peptide sequences.

    Parameters
    ----------
    peptide1 : List[str]
        The first tokenized peptide sequence to be compared.
    peptide2 : List[str]
        The second tokenized peptide sequence to be compared.
    aa_dict : Dict[str, float]
        Mapping of amino acid tokens to their mass values.
    cum_mass_threshold : float
        Mass threshold in Dalton to accept cumulative mass-matching amino acid
        sequences.
    ind_mass_threshold : float
        Mass threshold in Dalton to accept individual mass-matching amino acids.

    Returns
    -------
    aa_matches : np.ndarray of length max(len(peptide1), len(peptide2))
        Boolean flag indicating whether each paired-up amino acid matches across
        both peptide sequences.
    pep_match : bool
        Boolean flag to indicate whether the two peptide sequences fully match.
    per_seq_aa_matches : Tuple[np.ndarray]
        TODO.
    """
    aa_matches = np.zeros(max(len(peptide1), len(peptide2)), np.bool_)
    aa_matches_1 = np.zeros(len(peptide1), np.bool_)
    aa_matches_2 = np.zeros(len(peptide2), np.bool_)
    iso_errs = []
    tols_aa = []
    tols_prefix = []

    # Find longest mass-matching prefix.
    i1, i2, cum_mass1, cum_mass2 = 0, 0, 0.0, 0.0
    while i1 < len(peptide1) and i2 < len(peptide2):
        aa_mass1 = get_token_mass(peptide1[i1])
        aa_mass2 = get_token_mass(peptide2[i2])
        tol_prefix = abs(mass_diff(cum_mass1 + aa_mass1, cum_mass2 + aa_mass2, False))
        tols_prefix.append(tol_prefix)
        tol = abs(mass_diff(aa_mass1, aa_mass2, False))
        tols_aa.append(tol)
        if (tol_prefix < cum_mass_threshold):
            match = (
                tol < ind_mass_threshold
            )
            if match:
                iso_err = get_isobaric_error(peptide1[i1], peptide2[i2])
                if iso_err:
                    iso_errs.append(iso_err)
    
            aa_matches[max(i1, i2)] = match
            aa_matches_1[i1] = match
            aa_matches_2[i2] = match

            i1, i2 = i1 + 1, i2 + 1
            cum_mass1, cum_mass2 = cum_mass1 + aa_mass1, cum_mass2 + aa_mass2

        elif cum_mass2 + aa_mass2 > cum_mass1 + aa_mass1:
            i1, cum_mass1 = i1 + 1, cum_mass1 + aa_mass1
        else:
            i2, cum_mass2 = i2 + 1, cum_mass2 + aa_mass2
    return aa_matches, aa_matches.all(), (aa_matches_1, aa_matches_2), iso_errs, [tols_aa, tols_prefix]


def aa_match_prefix_suffix(
    peptide1: List[str],
    peptide2: List[str],
    aa_dict: Dict[str, float],
    cum_mass_threshold: float = 50,
    ind_mass_threshold: float = 20,
) -> Tuple[np.ndarray, bool, Tuple[np.ndarray]]:
    """
    Find the matching prefix and suffix amino acids between two peptide
    sequences.

    Parameters
    ----------
    peptide1 : List[str]
        The first tokenized peptide sequence to be compared.
    peptide2 : List[str]
        The second tokenized peptide sequence to be compared.
    aa_dict : Dict[str, float]
        Mapping of amino acid tokens to their mass values.
    cum_mass_threshold : float
        Mass threshold in Dalton to accept cumulative mass-matching amino acid
        sequences.
    ind_mass_threshold : float
        Mass threshold in Dalton to accept individual mass-matching amino acids.

    Returns
    -------
    aa_matches : np.ndarray of length max(len(peptide1), len(peptide2))
        Boolean flag indicating whether each paired-up amino acid matches across
        both peptide sequences.
    pep_match : bool
        Boolean flag to indicate whether the two peptide sequences fully match.
    per_seq_aa_matches : Tuple[np.ndarray]
        TODO.
    """
    # Find longest mass-matching prefix.
    aa_matches, pep_match, (aa_matches_1, aa_matches_2), iso_errs, tols = aa_match_prefix(
        peptide1, peptide2, aa_dict, cum_mass_threshold, ind_mass_threshold
    )

    # No need to evaluate the suffixes if the sequences already fully match.
    if pep_match:
        return aa_matches, pep_match, (aa_matches_1, aa_matches_2), iso_errs, tols

    # Find longest mass-matching suffix.
    i1, i2 = len(peptide1) - 1, len(peptide2) - 1
    i_stop = np.argwhere(~aa_matches)[0]
    cum_mass1, cum_mass2 = 0.0, 0.0 
    iso_errs = []
    tols_aa = []
    tols_suffix = []
    
    while i1 >= i_stop and i2 >= i_stop:
        aa_mass1 = get_token_mass(peptide1[i1])
        aa_mass2 = get_token_mass(peptide2[i2])
        tol_suffix = abs(mass_diff(cum_mass1 + aa_mass1, cum_mass2 + aa_mass2, False))
        tols_suffix.append(tol_suffix)
        tol_aa = abs(mass_diff(aa_mass1, aa_mass2, False))
        tols_aa.append(tol_aa)
        if (tol_suffix < cum_mass_threshold):

            match = tol_aa < ind_mass_threshold
            if match:
                iso_err = get_isobaric_error(peptide1[i1], peptide2[i2])
                if iso_err:
                    iso_errs.append(iso_err)
            
            aa_matches[max(i1, i2)] = match
            aa_matches_1[i1] = match
            aa_matches_2[i2] = match

            i1, i2 = i1 - 1, i2 - 1
            cum_mass1, cum_mass2 = cum_mass1 + aa_mass1, cum_mass2 + aa_mass2

        elif cum_mass2 + aa_mass2 > cum_mass1 + aa_mass1:
            i1, cum_mass1 = i1 - 1, cum_mass1 + aa_mass1
        else:
            i2, cum_mass2 = i2 - 1, cum_mass2 + aa_mass2

    # tols_aa = list(reversed(tols_aa))
    # tols_suffix = list(reversed(tols_aa))

    tols[0] += tols_aa
    tols[1] += tols_suffix

    return aa_matches, aa_matches.all(), (aa_matches_1, aa_matches_2), iso_errs, [tols_aa, tols_suffix]


def get_isobaric_error(aa1, aa2):
    '''Isobaric error type, boolean (True if resolved)'''
    iso_err_dict = {
        # Deemed irresolvable except for w-ions
        "I": "L",
        "L": "I",

        # semi-isobaric (slight diff in mass)
        "M": "F",
        "F": "M",
        "K": "Q",
        "Q": "K",

        # Sample prep isobarics
        "D": "N",
        "N": "D",
        "E": "Q",
        "Q": "E",
    }
    for i1, i2 in iso_err_dict.items():
        if aa1[0]==i1 and aa2[0]==i2:
            return (i1+"/"+i2, False)
        elif aa1[0]==i1 and aa2[0]==i1:
            return (i1+"/"+i2, True)


def aa_match(
    peptide1: List[str],
    peptide2: List[str],
    aa_dict: Dict[str, float] = {},
    cum_mass_threshold: float = 50,
    ind_mass_threshold: float = 20,
    mode: str = "best",
) -> Tuple[np.ndarray, bool, Tuple[np.ndarray]]:
    """
    Find the matching amino acids between two peptide sequences.

    Parameters
    ----------
    peptide1 : List[str]
        The first tokenized peptide sequence to be compared.
    peptide2 : List[str]
        The second tokenized peptide sequence to be compared.
    aa_dict : Dict[str, float]
        Mapping of amino acid tokens to their mass values.
    cum_mass_threshold : float
        Mass threshold in Dalton to accept cumulative mass-matching amino acid
        sequences.
    ind_mass_threshold : float
        Mass threshold in Dalton to accept individual mass-matching amino acids.
    mode : {"best", "forward", "backward"}
        The direction in which to find matching amino acids.

    Returns
    -------
    aa_matches : np.ndarray of length max(len(peptide1), len(peptide2))
        Boolean flag indicating whether each paired-up amino acid matches across
        both peptide sequences.
    pep_match : bool
        Boolean flag to indicate whether the two peptide sequences fully match.
    per_seq_aa_matches : Tuple[np.ndarray]
        TODO.
    """
    if mode == "best":
        return aa_match_prefix_suffix(
            peptide1, peptide2, aa_dict, cum_mass_threshold, ind_mass_threshold
        )
    elif mode == "forward":
        return aa_match_prefix(
            peptide1, peptide2, aa_dict, cum_mass_threshold, ind_mass_threshold
        )
    elif mode == "backward":
        aa_matches, pep_match, (aa_matches_1, aa_matches_2), iso_errs, tols = aa_match_prefix(
            list(reversed(peptide1)),
            list(reversed(peptide2)),
            aa_dict,
            cum_mass_threshold,
            ind_mass_threshold,
        )
        return (
            aa_matches[::-1],
            pep_match,
            (aa_matches_1[::-1], aa_matches_2[::-1]),
            iso_errs,
            tols
        )
    else:
        raise ValueError("Unknown evaluation mode")


def aa_match_batch(
    peptides1: Iterable,
    peptides2: Iterable,
    aa_dict: Dict[str, float],
    cum_mass_threshold: float = 50,
    ind_mass_threshold: float = 20,
    mode: str = "best",
) -> Tuple[List[Tuple[np.ndarray, bool]], int, int]:
    """
    Find the matching amino acids between multiple pairs of peptide sequences.

    Parameters
    ----------
    peptides1 : Iterable
        The first list of peptide sequences to be compared.
    peptides2 : Iterable
        The second list of peptide sequences to be compared.
    aa_dict : Dict[str, float]
        Mapping of amino acid tokens to their mass values.
    cum_mass_threshold : float
        Mass threshold in Dalton to accept cumulative mass-matching amino acid
        sequences.
    ind_mass_threshold : float
        Mass threshold in Dalton to accept individual mass-matching amino acids.
    mode : {"best", "forward", "backward"}
        The direction in which to find matching amino acids.

    Returns
    -------
    aa_matches_batch : List[Tuple[np.ndarray, bool, Tuple[np.ndarray]]]
        For each pair of peptide sequences: (i) boolean flags indicating whether
        each paired-up amino acid matches across both peptide sequences, (ii)
        boolean flag to indicate whether the two peptide sequences fully match,
        (iii) TODO.
    n_aa1: int
        Total number of amino acids in the first list of peptide sequences.
    n_aa2: int
        Total number of amino acids in the second list of peptide sequences.
    """
    aa_matches_batch, n_aa1, n_aa2 = [], 0, 0
    for peptide1, peptide2 in zip(peptides1, peptides2):
        # Split peptides into individual AAs if necessary.
        if isinstance(peptide1, str):
            peptide1 = re.split(r"(?<=.)(?=[A-Z])", peptide1)
        if isinstance(peptide2, str):
            peptide2 = re.split(r"(?<=.)(?=[A-Z])", peptide2)
        n_aa1, n_aa2 = n_aa1 + len(peptide1), n_aa2 + len(peptide2)
        aa_matches_batch.append(
            aa_match(
                peptide1,
                peptide2,
                aa_dict,
                cum_mass_threshold,
                ind_mass_threshold,
                mode,
            )
        )
    return aa_matches_batch, n_aa1, n_aa2


def aa_match_metrics(
    aa_matches_batch: List[Tuple[np.ndarray, bool]],
    n_aa_true: int,
    n_aa_pred: int,
) -> Tuple[float, float, float]:
    """
    Calculate amino acid and peptide-level evaluation metrics.

    Parameters
    ----------
    aa_matches_batch : List[Tuple[np.ndarray, bool]]
        For each pair of peptide sequences: (i) boolean flags indicating whether
        each paired-up amino acid matches across both peptide sequences, (ii)
        boolean flag to indicate whether the two peptide sequences fully match.
    n_aa_true: int
        Total number of amino acids in the true peptide sequences.
    n_aa_pred: int
        Total number of amino acids in the predicted peptide sequences.

    Returns
    -------
    aa_precision: float
        The number of correct AA predictions divided by the number of predicted
        AAs.
    aa_recall: float
        The number of correct AA predictions divided by the number of true AAs.
    pep_precision: float
        The number of correct peptide predictions divided by the number of
        peptides.
    """
    n_aa_correct = sum(
        [aa_matches[0].sum() for aa_matches in aa_matches_batch]
    )
    aa_precision = n_aa_correct / (n_aa_pred + 1e-8)
    aa_recall = n_aa_correct / (n_aa_true + 1e-8)
    pep_precision = sum([aa_matches[1] for aa_matches in aa_matches_batch]) / (
        len(aa_matches_batch) + 1e-8
    )
    return float(aa_precision), float(aa_recall), float(pep_precision)


def aa_precision_recall(
    aa_scores_correct: List[float],
    aa_scores_all: List[float],
    n_aa_total: int,
    threshold: float,
) -> Tuple[float, float]:
    """
    Calculate amino acid level precision and recall at a given score threshold.

    Parameters
    ----------
    aa_scores_correct : List[float]
        Amino acids scores for the correct amino acids predictions.
    aa_scores_all : List[float]
        Amino acid scores for all amino acids predictions.
    n_aa_total : int
        The total number of amino acids in the predicted peptide sequences.
    threshold : float
        The amino acid score threshold.

    Returns
    -------
    aa_precision: float
        The number of correct amino acid predictions divided by the number of
        predicted amino acids.
    aa_recall: float
        The number of correct amino acid predictions divided by the total number
        of amino acids.
    """
    n_aa_correct = sum([score > threshold for score in aa_scores_correct])
    n_aa_predicted = sum([score > threshold for score in aa_scores_all])
    return n_aa_correct / n_aa_predicted, n_aa_correct / n_aa_total


def convert_peptidoform(peptidoform):
    out = []
    n_mod = peptidoform.properties["n_term"]
    if n_mod is None:
        n_mod = [None]

    # If there is an N-terminal mod, this is seperately tokenized.
    else:
        out.append(("", n_mod))

    for i, aa_mod in enumerate(peptidoform):
        aa, mod = aa_mod
        if mod is None:
            mod = [mod]

        out.append((aa, mod))
    return out

def calculate_prc(scores_correct, scores_all, total_predicted, threshold):

    c = sum([score > threshold for score in scores_correct])
    ci = sum([score > threshold for score in scores_all])
    u = total_predicted-ci

    # precision
    precision = c/ci

    # recall (This is an alternative definition and the line will stop at x=y)
    recall = c/total_predicted

    # coverage
    coverage = ci/total_predicted

    return precision, recall, coverage



def get_prc_curve(t, total_predicted):

    prs = []
    recs = []
    covs = []

    for threshold in np.linspace(t.score.max(), t.score.min(), 1000):
        if np.isnan(threshold):
            continue
        
        pr, rec, cov = calculate_prc(
            scores_correct=t[t.match].score.to_numpy(),
            scores_all=t.score.to_numpy(),
            total_predicted=total_predicted,
            threshold=threshold
        )
        prs.append(pr)
        recs.append(rec)
        covs.append(cov)

    return pd.DataFrame(
        {'precision': prs,
         'recall': recs,
         'coverage': covs}
    )

def annotate_stacked_histplot(ax):
    # Get the bars' heights
    patches = ax.patches

    # Grouping bars by bins
    heights_by_bin = {}

    # Loop through each patch and store heights by bin
    for patch in patches:
        # The bin location is the x-coordinate of the patch
        bin_location = patch.get_x()
        
        # Store height based on bin location
        if bin_location not in heights_by_bin:
            heights_by_bin[bin_location] = []
        heights_by_bin[bin_location].append(patch.get_height())

    # Annotate percentages
    for bin_location, heights in heights_by_bin.items():
        total_height = sum(heights)
        
        # Calculate percentage and annotate each stack if greater than 5%
        cumulative_height = 0
        for height in heights:
            percentage = (height / total_height) * 100
            
            if percentage > 5:
                # Get the center x-coordinate of the patch
                x = bin_location + 0.25  # Adjust based on binwidth
                
                # Annotate the center of each bar in the stack
                y = cumulative_height + height / 2  # Midpoint of the bar
                
                # Add annotation to the plot
                ax.text(x, y, f'{percentage:.0f}%', ha='center', va='center', color='black', fontsize=8)
            
            # Update cumulative height for stacking annotations
            cumulative_height += height

def get_refinement_error_tables(run, base_engine, score_metadata='ms2rescore', post_processors=["Spectralis", "InstaNovo+"]):
    error_types = {
        post_processor: {"Unchanged": {
            "match": 0,
            "isobaric_aa": 0,
            "isobaric_peak": 0,
            "Higher score": 0,
            "Lower score": 0
        },
        "Refined_base": {
            "match": 0,
            "isobaric_aa": 0,
            "isobaric_peak": 0,
            "Higher score": 0,
            "Lower score": 0
        },
        "Refined": {
            "match": 0,
            "isobaric_aa": 0,
            "isobaric_peak": 0,
            "Higher score": 0,
            "Lower score": 0
        }} for post_processor in post_processors
    }

    for spectra in run.spectra.values():
        psms = spectra.get_psms_by_engine(base_engine)
        if len(psms) == 0:
            continue
        psm = psms[0]
        for r_name, (r_psm, unchanged) in psm.refinement.items():
            
            if unchanged:
                error_type = psm.evaluation[score_metadata].error_type
                error_types[r_name]["Unchanged"][error_type] += 1
            else:
                error_type = r_psm.evaluation[score_metadata].error_type
                error_types[r_name]["Refined"][error_type] += 1

                error_type = psm.evaluation[score_metadata].error_type
                error_types[r_name]["Refined_base"][error_type] += 1

    error_table = {
        post_processor: pd.DataFrame(error_types[post_processor]) for post_processor in post_processors
    }
    return error_table

def plot_refinement_error_table(error_table, fig, ax):

    refined_pct = {}
    for i, (post_processor, t) in enumerate(error_table.items()):
        refined_pct[post_processor] = t.sum()["Refined"] / t.sum()[["Unchanged", "Refined"]].sum()

        t.T.plot(
            kind="bar",
            stacked=True,
            ax=ax[i],
            title="Post-processing with {}\nChanged {:.2f}% of PSMs".format(post_processor, refined_pct[post_processor]*100)
        )

        annotate_stacked_histplot(ax[i])
        ax[i].tick_params(axis='x', rotation=0)

def get_match_score_table(run, engine, score_metadata, eval_score_metadata="ms2rescore", refiner=None, return_table=False):

    pr_table = {
        "spectrum_id": [],
        "score": [],
        "match": []
    }
    missing = 0

    for spectrum in run.spectra.values():
        
        psm = spectrum.get_psms_by_engine(engine)[0]

        if refiner is not None:
            if refiner not in psm.refinement.keys():
                missing+=1
                continue
            if psm.refinement[refiner][1]:
                pr_table["match"].append(psm.evaluation[eval_score_metadata].error_type=="match")
                pr_table["score"].append(psm.scores.get_score(score_metadata))
            else:
                psm_refinement = psm.refinement[refiner][0]
                pr_table["match"].append(psm_refinement.evaluation[eval_score_metadata].error_type=="match")
                pr_table["score"].append(psm_refinement.scores.get_score(score_metadata))
            pr_table["spectrum_id"].append(spectrum.spectrum_id)

        else:
            pr_table["match"].append(psm.evaluation[eval_score_metadata].error_type=="match")
            pr_table["score"].append(psm.scores.get_score(score_metadata))
            pr_table["spectrum_id"].append(spectrum.spectrum_id)
    pr_table = pd.DataFrame(pr_table)

    if return_table:
        return pr_table

    return precision_recall_curve(
        pr_table.match.tolist(),
        pr_table.score.tolist()), missing


def plot_precision_recall_refinement(
        denovo_seqs,
        spectralis_seqs,
        instanovo_seqs,
        denovo_name,
        ax,
        max_points_per_plot=1000
):

    # Sample color mappings
    denovo_color = "r"
    spectralis_color = "b"
    instanovo_color = "g"

    # Line style mappings
    score_type_styles = {
        denovo_name: "solid",
        "InstaNovo+": "solid",
        "Spectralis": "dashdot",
        "ms2rescore": "dotted"
    }

    # Create figure and axis

    print("denovo")
    # Loop over the denovo sequences
    for score_type, pr_dict in denovo_seqs.items():
        print("plotting", score_type)
        sns.lineplot(
            x=pr_dict["recall"],
            y=pr_dict["precision"],
            label=score_type,   # Label for style legend
            linestyle=score_type_styles[score_type],  # Set line style
            color=denovo_color,  # Set color
            alpha=.5,
            ax=ax
        )

    print("spectralis")
    # Loop over the spectralis sequences
    for score_type, pr_dict in spectralis_seqs.items():
        print("plotting", score_type)
        sns.lineplot(
            x=pr_dict["recall"],
            y=pr_dict["precision"],
            label=score_type,   # Label for style legend
            linestyle=score_type_styles[score_type],  # Set line style
            color=spectralis_color,  # Set color
            alpha=.5,
            ax=ax
        )

    print("instanovo+")
    # Loop over the instanovo sequences
    for score_type, pr_dict in instanovo_seqs.items():
        print("plotting", score_type)
        sns.lineplot(
            x=pr_dict["recall"],
            y=pr_dict["precision"],
            label=score_type,   # Label for style legend
            linestyle=score_type_styles[score_type],  # Set line style
            color=instanovo_color,  # Set color
            alpha=.5,
            ax=ax
        )

    # Create separate legends
    # First, get handles and labels from the plot
    handles, labels = ax.get_legend_handles_labels()

    # Separate out color-based and style-based labels for different legends
    color_legend_labels = [denovo_name, 'Spectralis', 'InstaNovo']
    style_legend_labels = [denovo_name, 'InstaNovo+', 'Spectralis', 'ms2rescore']

    # Add legends outside the plot
    color_handles = [plt.Line2D([0], [0], color=c, lw=3) for c in [denovo_color, spectralis_color, instanovo_color]]
    style_handles = [plt.Line2D([0], [0], linestyle=score_type_styles[label], color='k', lw=1) for label in style_legend_labels]
    leg1 = ax.legend(color_handles, color_legend_labels, title="Sequence", loc='upper left', bbox_to_anchor=(1, 1))  # Outside top-right
    leg2 = ax.legend(style_handles, style_legend_labels, title="Score Type", loc='lower left', bbox_to_anchor=(1, 0))  # Outside bottom-right

    ax.add_artist(leg1)
    ax.add_artist(leg2)
    ax.legend(color_handles, color_legend_labels, title="Dataset", loc='upper left', bbox_to_anchor=(1, 1))  # Outside top-right

    # Adjust plot to make space for legends
    plt.tight_layout(rect=[0, 0, 0.85, 1])  # Adjust space for legends outside the plot

    # Ensure both legends are shown
    plt.show()


def load_seq_score_dicts(
        run,
        denovo_name,
        total_predicted
):
    denovo_seqs = {
        denovo_name: {},
        "ms2rescore": {},
        "Spectralis": {}
    }

    for score_metadata in denovo_seqs.keys():
        t = get_match_score_table(
            run, engine=denovo_name, score_metadata=score_metadata, return_table=True
        )       
        prc = get_prc_curve(t, total_predicted)

        denovo_seqs[score_metadata] = {
            'precision': prc.precision.to_numpy(),
            'recall': prc.recall.to_numpy(),
            'coverage': prc.coverage.to_numpy()
        }

    spectralis_seqs = {
        "Spectralis": [],
        "ms2rescore": []
    }

    for score_metadata in spectralis_seqs.keys():
        t = get_match_score_table(
            run, engine=denovo_name, score_metadata=score_metadata, refiner="Spectralis", return_table=True
        )
        
        prc = get_prc_curve(t, total_predicted)
        spectralis_seqs[score_metadata] = {
            'precision': prc.precision.to_numpy(),
            'recall': prc.recall.to_numpy(),
            'coverage': prc.coverage.to_numpy()
        }
    
    instanovo_seqs = {
        "InstaNovo+": [],
        "ms2rescore": []
    }

    for score_metadata in instanovo_seqs.keys():
        t = get_match_score_table(
            run, engine=denovo_name, score_metadata=score_metadata, refiner="InstaNovo+", return_table=True
        )
        prc = get_prc_curve(t, total_predicted)

        instanovo_seqs[score_metadata] = {
            'precision': prc.precision.to_numpy(),
            'recall': prc.recall.to_numpy(),
            'coverage': prc.coverage.to_numpy()
        }

    seq_score_dict = {
        "denovo_seqs": denovo_seqs,
        "spectralis_seqs": spectralis_seqs,
        "instanovo_seqs": instanovo_seqs
    }
    return seq_score_dict