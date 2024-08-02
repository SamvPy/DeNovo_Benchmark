import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pyteomics import mgf
import spectrum_utils.plot as sup
import spectrum_utils.spectrum as sus

def plot_metrics(metrics_dict, title: str):
    barplot = sns.barplot(
        pd.DataFrame(metrics_dict).reset_index().melt(id_vars="index"),
        y="value",
        x="index",
        hue="variable"
    )

    for p in barplot.patches:
        height = p.get_height()
        barplot.text(
            p.get_x() + p.get_width() / 2.0,  # X-coordinate: center of the bar
            height,  # Y-coordinate: height of the bar
            f'{height*100:.1f}%',  # Text to be displayed
            ha='center',  # Horizontal alignment
            va='bottom'   # Vertical alignment
        )

    # Example title: "% spectra extra predicted by de novo tools above threshold (FILTERED)"
    plt.title()
    plt.ylim((0,1.1))


def plot_spectrum(mgf_path, spectrum_id, peptide, fragment_tol_mass=50, fragment_tol_mode="ppm", ion_types="by", plot=True):
    mgf_file = mgf.read(mgf_path)
    spectrum = mgf_file.get_by_id(spectrum_id)

    spectrum_su = sus.MsmsSpectrum(
        identifier=spectrum["params"]["title"],
        precursor_mz=spectrum["params"]["pepmass"][0],
        precursor_charge=spectrum["params"]["charge"][0],
        mz=spectrum["m/z array"],
        intensity=spectrum["intensity array"],
        retention_time=spectrum["params"]["rtinseconds"]
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
