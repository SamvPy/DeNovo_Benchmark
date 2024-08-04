"""Functions for visualizing metrics, and mass spectra."""

from typing import Literal

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import spectrum_utils.plot as sup
import spectrum_utils.spectrum as sus
from pyteomics import mgf


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
    mgf_path: str,
    spectrum_id: str,
    peptide: str,
    fragment_tol_mass: float = 50.0,
    fragment_tol_mode: Literal["ppm", "Da"] = "ppm",
    ion_types: str = "by",
    plot: bool = True,
) -> sus.MsmsSpectrum:
    """
    Plot a mass spectrum with peak annotations.

    Parameters
    ----------
    mgf_path: str
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
    mgf_file = mgf.read(mgf_path)
    spectrum = mgf_file.get_by_id(spectrum_id)

    spectrum_su = sus.MsmsSpectrum(
        identifier=spectrum["params"]["title"],
        precursor_mz=spectrum["params"]["pepmass"][0],
        precursor_charge=spectrum["params"]["charge"][0],
        mz=spectrum["m/z array"],
        intensity=spectrum["intensity array"],
        retention_time=spectrum["params"]["rtinseconds"],
    )
    spectrum_su = spectrum_su.annotate_proforma(
        proforma_str=peptide,
        fragment_tol_mass=fragment_tol_mass,
        fragment_tol_mode=fragment_tol_mode,
        ion_types=ion_types,
    )

    if plot:
        sup.spectrum(spectrum_su)
    return spectrum_su
