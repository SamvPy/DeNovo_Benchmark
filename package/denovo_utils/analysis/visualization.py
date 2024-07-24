import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

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