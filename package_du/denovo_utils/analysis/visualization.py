"""Functions for visualizing metrics, and mass spectra."""

from typing import Literal

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import spectrum_utils.plot as sup
import spectrum_utils.spectrum as sus
import numpy as np
from scipy import stats
from pyteomics import mgf
from pyteomics.mzml import PreIndexedMzML


def plot_metrics(metrics_dict: dict, title: str) -> None:
    """
    Plot metrics in a barplot.

    Parameters
    ----------
    metrics_dict: dict
        A dictionary containing the keys as metric label, and values the scores.
    title: str
        Title to give to the barplot
    """
    barplot = sns.barplot(
        pd.DataFrame(metrics_dict).reset_index().melt(id_vars="index"),
        y="value",
        x="index",
        hue="variable",
    )

    for p in barplot.patches:
        height = p.get_height()
        barplot.text(
            p.get_x() + p.get_width() / 2.0,  # X-coordinate: center of the bar
            height,  # Y-coordinate: height of the bar
            f"{height*100:.1f}%",  # Text to be displayed
            ha="center",  # Horizontal alignment
            va="bottom",  # Vertical alignment
        )

    # Example title: "% spectra extra predicted by de novo tools
    # above threshold (FILTERED)"
    plt.title(title)
    plt.ylim((0, 1.1))


def plot_spectrum(
    peak_path: str,
    spectrum_id: str,
    peptide: str,
    peak_format: str = 'mzml',
    fragment_tol_mass: float = 20,
    fragment_tol_mode: Literal["ppm", "Da"] = "ppm",
    ion_types: str = "by",
    max_ion_charge: int = 2,
    plot: bool = True,
) -> sus.MsmsSpectrum:
    """
    Plot a mass spectrum with peak annotations.

    Parameters
    ----------
    peak_path: str
        Path to the mgf file containing the spectrum to plot.
    spectrum_id: str
        The spectrum_id of the spectrum in the mgf file.
    peptide: str (in proforma format)
        Used for peak annotations.
    fragment_tol_mass: float (default=50)
        Fragment tolerance for peak annotation.
    fragment_tol_mode: Literal['Da', 'ppm] (default=ppm)
        Mode of fragment mass tolerance, either in Dalton or parts-per-million.
    plot: bool (default=True)
        Whether to plot the annotated spectrum.

    Return
    ------
    MsmsSpectrum
        The annotated mass spectrum in spectrum_utils format.
    """
    if peak_format=='mzml':
        mzml_file = PreIndexedMzML(peak_path)
        spectrum = mzml_file.get_by_id(spectrum_id)

        title = spectrum['id']
        selected_ion = spectrum['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]
        precursor_mz = float(selected_ion['selected ion m/z'])
        retention_time = spectrum['scanList']['scan'][0]['scan start time'] * 60
        precursor_charge = int(selected_ion['charge state'])


    elif peak_format=='mgf':
        mgf_file = mgf.read(peak_path)
        spectrum = mgf_file.get_by_id(spectrum_id)

        title = spectrum["params"]["title"]
        precursor_mz = spectrum["params"]["pepmass"][0]
        precursor_charge = spectrum["params"]["charge"][0]
        retention_time = spectrum["params"]["rtinseconds"]

    else:
        raise Exception('Peak format {} not supported. Only mzml or mgf.'.format(peak_format))

    spectrum_su = sus.MsmsSpectrum(
        identifier=title,
        precursor_mz=precursor_mz,
        precursor_charge=precursor_charge,
        mz=spectrum["m/z array"],
        intensity=spectrum["intensity array"],
        retention_time=retention_time,
    )

    spectrum_su = spectrum_su.annotate_proforma(
        proforma_str=peptide,
        fragment_tol_mass=fragment_tol_mass,
        fragment_tol_mode=fragment_tol_mode,
        ion_types=ion_types,
        max_ion_charge=max_ion_charge
    )

    if plot:
        sup.spectrum(spectrum_su)
    return spectrum_su


def plot_gmm_fit(scores: list, gmm_dict: dict, threshold: float):

    x = np.linspace(-3, 0, 5000)
    scores_arr = np.array(scores)

    # Calculate the PDF of each normal distribution
    pdf1 = gmm_dict["weight"][0] * stats.norm.pdf(x, gmm_dict["mean"][0], gmm_dict["std"][0])
    pdf2 = gmm_dict["weight"][1] * stats.norm.pdf(x, gmm_dict["mean"][1], gmm_dict["std"][1])

    kde = stats.gaussian_kde(scores_arr)
    kde_obs = kde(x)

    # Plot the individual distributions
    plt.plot(x, pdf1, label='Pseudo decoys', linestyle=":")
    plt.plot(x, pdf2, label='Pseudo targets', linestyle=":")

    # Plot the mixture distribution
    plt.plot(x, pdf1 + pdf2, label='Weighted Mixture', linestyle='--')
    plt.plot(x, kde_obs, label="Original data", linestyle="solid", color="grey")

    plt.axvline(threshold, alpha=.75)

    plt.legend()
    plt.title('GMM fit with threshold at {:.2f}'.format(threshold))
    plt.show()


def plot_score_distribution(run, spectrum_id: str, score_name='score_ms2rescore_features'):
    spectrum = run.get_spectrum(spectrum_id)
    
    ms2rescore_gt = spectrum.psm_gt.scores.get_score(score_name)
    ms2rescore_candidates = [
        psm.scores.get_score(score_name) for psm in
        spectrum.psm_candidates if psm.engine_name in ['Terminal-variant', 'Ambiguous']
    ]

    candidate_type = [
        psm.engine_name for psm in spectrum.psm_candidates if psm.engine_name in [
            'Terminal-variant', 'Ambiguous'
        ]
    ]

    engine_output = [
        psm.scores.get_score(score_name) for psm in 
        spectrum.psm_candidates if psm.engine_name not in ['Terminal-variant', 'Ambiguous']
    ]

    sns.histplot(
        x=ms2rescore_candidates,
        hue=candidate_type,
        binwidth=.2,
        kde=True
    )

    for engine_score in engine_output:
        plt.axvline(engine_score, c='g')
    plt.axvline(ms2rescore_gt, c='r')

    plt.xlim(
        min([-10]+[ms2rescore_gt]+ms2rescore_candidates),
        max([10]+[ms2rescore_gt]+ms2rescore_candidates)
    )