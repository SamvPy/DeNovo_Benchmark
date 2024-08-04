"""
Functionalities for Gaussian Mixture Modelling of score distributions.

Rewrite this part!
"""

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats
from sklearn.mixture import GaussianMixture as GMM


def get_ggm_n_clusters(scores: list, search_engine: str = "") -> None:
    """
    Extract the optimal number of clusters (evaluate up to 10 components).

    Parameters
    ----------
    scores: list
        Scores for which to fit a GMM.
    search_engine: str
        The label of the search engine. Only used for aesthetics of plot.
    """
    # import libraries (some are for cosmetics)

    mpl.rcParams["axes.linewidth"] = 1.5
    mpl.rcParams.update(
        {"font.size": 15, "font.family": "STIXGeneral", "mathtext.fontset": "stix"}
    )

    # create the data as in @Meng's answer
    x = scores
    x = x.reshape(-1, 1)

    # first of all, let's confirm the optimal number of components
    bics = []
    min_bic = 0
    counter = 1
    for i in range(10):  # test the AIC/BIC metric between 1 and 10 components
        gmm = GMM(
            n_components=counter, max_iter=1000, random_state=0, covariance_type="full"
        )
        _ = gmm.fit(x).predict(x)
        bic = gmm.bic(x)
        bics.append(bic)
        if bic < min_bic or min_bic == 0:
            min_bic = bic
            opt_bic = counter
            counter = counter + 1

    # plot the evolution of BIC/AIC with the number of components
    fig = plt.figure(figsize=(10, 4))
    _ = fig.add_subplot(1, 2, 1)
    # Plot 1
    plt.plot(np.arange(1, 11), bics, "o-", lw=3, c="black", label="BIC")
    plt.legend(frameon=False, fontsize=15)
    plt.xlabel("Number of components", fontsize=20)
    plt.ylabel("Information criterion", fontsize=20)
    plt.xticks(np.arange(0, 11, 2))
    plt.title(f"{search_engine}\tOpt. components = {opt_bic}", fontsize=20)


def assign_ggm_clusters(
    scores: list, n_clusters: int = 4
) -> tuple[
    np.typing.ArrayLike, np.typing.ArrayLike, np.typing.ArrayLike, np.typing.ArrayLike
]:
    """
    Assign cluster labels by modeling score list with GMM with n_cluster components.

    Parameters
    ----------
    scores: list
        List of scores
    n_clusters: int (default=4)
        Number of components to model in GMM.

    Returns
    -------
    labels: np.array
        Label assignment to cluster based on max probability.
    mean: np.array
        Array of means for the clusters.
    covs:
        Covariates of the clusters.
    weights:
        Weights for the clusters.
    """
    # create GMM model object
    x = scores
    x = x.reshape(-1, 1)
    gmm = GMM(
        n_components=n_clusters, max_iter=1000, random_state=10, covariance_type="full"
    )

    # find useful parameters
    mean = gmm.fit(x).means_.flatten()
    covs = gmm.fit(x).covariances_.flatten()
    weights = gmm.fit(x).weights_.flatten()

    return gmm.predict(scores.reshape(-1, 1)), mean, covs, weights


def plot_cluster_psmtype(
    df: pd.DataFrame,
    mean: np.typing.ArrayLike,
    covs: np.typing.ArrayLike,
    weights: np.typing.ArrayLike,
    n_clusters: int = 4,
    score_col: str = "score",
    psm_type_col: str = "psm_type",
) -> None:
    """
    Plot a KDE-plot to see how the GMM created clusters from a score distribution.

    Parameters
    ----------
    df: pd.DataFrame
        Dataframe with columns score_col and psm_type_col.
    mean:
        Means of the clusters
    covs:
        Covariates of the clusters.
    weights:
        Weights of the clusters.
    n_clusters: int (default=4)
        The number of clusters.
        Number of figures plotted corresponds with this number.
    score_col: str (default=score)
        The column name of the scores to plot in the KDE-plot.
    psm_type_col: str (default=psm_type)
        The column name of psm-types (colors in plot).
    """
    _, ax = plt.subplots(1, n_clusters + 1, figsize=(25, 6))

    x_axis = np.linspace(-1, 1, 1000)
    n_gaussians = n_clusters

    colors = sns.color_palette("tab20")[:n_gaussians]
    for color, i in zip(colors, list(range(n_gaussians))):
        ax[0].plot(
            x_axis,
            weights[i] * stats.norm.pdf(x_axis, mean[i], np.sqrt(covs[i])),
            c=color,
        )
    _ = ax[0].hist(df[score_col], density=True, bins=1000)

    for cluster in range(0, n_clusters):
        subset = df[df.gmm_cluster == cluster].sort_values(psm_type_col)
        sns.kdeplot(subset, x=score_col, hue=psm_type_col, ax=ax[cluster + 1])
        ax[cluster + 1].set_title(
            "Cluster {}\nMean: {:.2f},  STD: {:.2f}".format(
                cluster, mean[cluster], np.sqrt(covs[cluster])
            )
        )
        ax[cluster + 1].set_xlim((df[score_col].min(), df[score_col].max()))


def filter_highest_cluster_psms(
    df: pd.DataFrame, source: str, score_col: str
) -> pd.DataFrame:
    """
    Filter a dataframe by extracting PSMs assigned to highest scoring cluster.

    Parameters
    ----------
    df: pd.DataFrame
        Dataframe to be filtered. Needs columns 'source', 'gmm_cluster', and score_col.
    source: str
        The label of the search engine to filter in 'source' column.
    score_col: str
        The column name in which scores are stored.
    Return
    ------
    pd.DataFrame
        Copied and filtered dataframe.
    """
    df = df.loc[df["source"] == source].copy()

    labels, _, _, _ = assign_ggm_clusters(df[score_col].to_numpy(), n_clusters=2)

    df["gmm_cluster"] = labels
    high_cluster_label = df.groupby("gmm_cluster")[score_col].max().idxmax()
    return df[df["gmm_cluster"] == high_cluster_label]
