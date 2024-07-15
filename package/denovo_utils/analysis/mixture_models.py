import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
from sklearn.mixture import GaussianMixture as GMM
import matplotlib as mpl
import seaborn as sns


def get_ggm_n_clusters(scores, search_engine=""):
    # import libraries (some are for cosmetics)

    mpl.rcParams['axes.linewidth'] = 1.5
    mpl.rcParams.update({'font.size': 15, 'font.family': 'STIXGeneral', 'mathtext.fontset': 'stix'})


    # create the data as in @Meng's answer
    x = scores
    x = x.reshape(-1, 1)

    # first of all, let's confirm the optimal number of components
    bics = []
    min_bic = 0
    counter=1
    for i in range (10): # test the AIC/BIC metric between 1 and 10 components
        gmm = GMM(n_components = counter, max_iter=1000, random_state=0, covariance_type = 'full')
        labels = gmm.fit(x).predict(x)
        bic = gmm.bic(x)
        bics.append(bic)
        if bic < min_bic or min_bic == 0:
            min_bic = bic
            opt_bic = counter
            counter = counter + 1


    # plot the evolution of BIC/AIC with the number of components
    fig = plt.figure(figsize=(10, 4))
    ax = fig.add_subplot(1,2,1)
    # Plot 1
    plt.plot(np.arange(1,11), bics, 'o-', lw=3, c='black', label='BIC')
    plt.legend(frameon=False, fontsize=15)
    plt.xlabel('Number of components', fontsize=20)
    plt.ylabel('Information criterion', fontsize=20)
    plt.xticks(np.arange(0,11, 2))
    plt.title(f'{search_engine}\tOpt. components = {opt_bic}', fontsize=20)


def assign_ggm_clusters(scores, n_clusters = 4):
    # create GMM model object
    x = scores
    x = x.reshape(-1, 1)
    gmm = GMM(n_components = n_clusters, max_iter=1000, random_state=10, covariance_type = 'full')

    # find useful parameters
    mean = gmm.fit(x).means_.flatten()
    covs  = gmm.fit(x).covariances_.flatten()
    weights = gmm.fit(x).weights_.flatten()

    print(
        f"Mean: {mean}\nCovariance: {covs}\nWeights: {weights}"
    )


    return gmm.predict(scores.reshape(-1, 1)), mean, covs, weights

def plot_cluster_psmtype(df, mean, covs, weights, n_clusters=4, score_col="score", psm_type_col="psm_type"):
    fig, ax = plt.subplots(1,n_clusters+1, figsize=(25,6))

    x_axis = np.linspace(-1, 1, 1000)
    n_gaussians = n_clusters


    colors = sns.color_palette("tab20")[:n_gaussians]
    for color, i in zip(colors, list(range(n_gaussians))):
        ax[0].plot(
            x_axis,
            weights[i]*stats.norm.pdf(x_axis, mean[i], np.sqrt(covs[i])),
            c=colors[i] 
        )
    p = ax[0].hist(df[score_col], density=True, bins=1000)

    for cluster in range(0,n_clusters):
        subset = df[df.gmm_cluster==cluster].sort_values(psm_type_col)
        sns.kdeplot(
            subset,
            x= score_col,
            hue=psm_type_col,
            ax=ax[cluster+1]
        )
        ax[cluster+1].set_title("Cluster {}\nMean: {:.2f},  STD: {:.2f}".format(cluster, mean[cluster], np.sqrt(covs[cluster])))
        ax[cluster+1].set_xlim((df[score_col].min(), df[score_col].max()))