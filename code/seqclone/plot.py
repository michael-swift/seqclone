from matplotlib import pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from scipy.stats import t


def plot_waterfall(df, adata_obs, gene, label):
    """ Strip plot where x axis is rank ordered by mean gene expression within a clone"""
    labels_testing = adata_obs
    fig, ax = plt.subplots(1, 1)
    labels_testing.loc[:, gene] = df[gene]
    order = labels_testing.groupby(label)[gene].mean().sort_values(ascending=False).index
    g = sns.stripplot(ax=ax, data=labels_testing, x=labels_testing[label], y=gene, order=order, color=None)
    # save_figure(fig, '{}_{}'.format(label, str(gene)))
    return g

def plotWaterfall2(gene, label, data, df, remove_zeroes, shuffle, title, save):
    
    # Only Clones
    if label == 'CLONE':
        _selector = data.CLONE.value_counts() > 1
        _selector = _selector.index[_selector == True]
        data = data[data.CLONE.isin(_selector)]
    # add gene expression data
    data.loc[:, gene] = df[gene]
    data = data[data[label].isin((data[label].value_counts() > 1).index)]
    sns.set(rc={'figure.figsize':(8, 3)})

    data = pd.DataFrame([data[label].to_list(), data[gene].to_list()])
    data = data.T
    data.columns = [label, 'gene_expression']
    if remove_zeroes == True:
        _df = data.groupby('CLONE').mean() > 0
        selector = _df[_df.gene_expression == True].index
        data = data[data.CLONE.isin(selector)]
    fig, ax = plt.subplots(1,1)
    if shuffle == True:
        np.random.shuffle(data.gene_expression.values)

    order = data.groupby(label)['gene_expression'].mean().sort_values().index
    y = data['gene_expression']
    ax = sns.stripplot(data = data, x = label, y = 'gene_expression', order = order, size = 4.5, linewidth = 0.5, alpha = 0.6, palette='dark')
    ax.set_xlabel('Lineage')
    ax.set_ylabel('log$_{10}$ CPM')
    ax.axhline(y.mean(), linestyle = 'dotted', c = 'k')
    ax.axhline(y.mean()+ y.std(), linestyle = '-.', c = 'grey')
    ax.axhline(y.mean() - y.std(), linestyle = '-.', c = 'grey')
    ax.tick_params(axis = 'x',labelbottom=False)
    name = title + "_" + gene
    ax.set_title(title + ' ' + gene )

    fig = ax.get_figure()
    
    if save == True: 
        save_figure(fig, name)
    return

def plot_swarm(df, adata_obs, gene, label):
    """ Strip plot where x axis is rank ordered by mean gene expression within a clone"""
    labels_testing = adata_obs
    fig, ax = plt.subplots(1, 1)
    labels_testing.loc[:, gene] = df[gene]
    order = labels_testing.groupby(label)[gene].mean().sort_values(ascending=False).index
    g = sns.swarmplot(ax=ax, data=labels_testing, x=labels_testing[label], y=gene, order=order, color=None)
    # save_figure(fig, '{}_{}'.format(label, str(gene)))
    return g


def plot_ci(df, adata, num_shuffles, gene, label, alpha):
    # This is time consuming to do twice, like this because I really am only plotting a subset of hits
    labels_testing = pd.merge(adata.obs[label], df[gene], left_index=True, right_index=True)
    tested_gene = []
    statistics = []
    # get the ordered means of the true labeling
    observedlabel_mean = labels_testing.groupby(label)[gene].mean().sort_values(ascending=False)
    ci_df = pd.DataFrame(observedlabel_mean)
    # shuffling loop
    for i in range(num_shuffles):
        # create copy out of superstition
        labels_testing_copy = labels_testing.copy(deep=True)
        # shuffle labels
        labels_testing_copy[label] = np.random.permutation(labels_testing_copy[label].values)
        # calculate mean gene expression of the shuffled
        shuffled_means = labels_testing_copy.groupby(label)[gene].mean()
        # add the shuffled mean to the ci_df
        ci_df = pd.merge(ci_df, shuffled_means, left_index=True, right_index=True)
    # First column is the true label means
    true_ordered_means = ci_df.iloc[:, 0]
    # rest of the columns are the shuffled means
    shuffled_means = ci_df.iloc[:, 1:]
    # Using T distribution
    T_low = shuffled_means.apply(lambda row: t.interval(alpha, row.shape[0]-1, loc=row.median(), scale=row.sem())[0], axis=1)
    T_high = shuffled_means.apply(lambda row: t.interval(alpha, row.shape[0]-1, loc=row.median(), scale=row.sem())[1], axis=1)
    lower_quantile = shuffled_means.quantile(q=0.025, axis=1)
    upper_quantile = shuffled_means.quantile(q=0.975, axis=1)

    # merge data for plotting
    """
    data = pd.merge(true_ordered_means, shuffled_means, left_index = True, right_index= True)
    data.reset_index(inplace = True)
    data[label] = data[label].astype('str')
    fig, ax = plt.subplots()
    x = data[label]
    # Frequentist confidence intervals
    f_lowci = data['lower_quant']
    f_upci = data['upper_quant']
    true_data = true_ordered_means
    # T distribution Quantiles
    g_lowci = data['lower']
    g_upci = data['upper']
    ax.plot(x, true_data, label = 'True Data Order')
    ax.fill_between(x, f_lowci, f_upci, alpha = alpha, color = 'k', label = 'CI using real quantiles')
    ax.fill_between(x, g_lowci, g_upci, alpha = alpha, color = 'r', label = 'CI using T distribution')
    plt.xlabel(label)
    plt.ylabel(gene + ' \n mean expression (log CPM)')
    plt.xticks()
    ax.legend() """
    data = 'k'
    return data, shuffled_means


def plot_test_hist(df, adata_obs, num_shuffles, gene, label):
    """ Plots the histogram of bootstrapped mean variances amongst clones vs. the true mean variance amongst clones"""
    """ Could another way to do this test be a KS test on all of the null distributions of mean variances generated vs.
     the true distribution of mean variances"""
    labels_testing = adata_obs.copy()
    tested_gene = []
    statistics = []
    # Select the gene to test
    labels_testing[gene] = df.loc[:, gene]
    # get the mean variance of the true clonal groupings
    observedlabel_var = labels_testing.groupby(label)[gene].var()
    mean_observedlabel_variance = observedlabel_var.mean()
    # Generate a null distribution of mean variances by shuffling clonal groupings
    mean_shuffled_variances = []
    for i in range(num_shuffles):
        # create copy, honestly don't think this is necessary, but only hurts performance
        labels_testing_copy = labels_testing.copy(deep=True)
        # shuffle labels
        labels_testing_copy.loc[:,label] = np.random.permutation(labels_testing_copy[label].values)
        # Calculate variances of shuffled clonal groups
        shuffled_variances = labels_testing_copy.groupby(label)[gene].var()
        # Calculate mean of variances of shuffled clonal groups
        mean_shuffled_variances.append(shuffled_variances.mean())
    # Plot
    fig, ax = plt.subplots(1,1)
    data = pd.Series(mean_shuffled_variances)
    xmax = data.max() + 0.2
    bins = np.linspace(0, xmax, 100)
    plt.hist(data, bins=bins, color='peru', alpha=0.5)
    plt.hist(data, bins=bins, color='peru', histtype='step')
    plt.axvline(mean_observedlabel_variance, 0, 1, c='midnightblue', ls='--')
    plt.yscale('log')
    plt.xscale('linear')
    plt.xlabel('Average ' + label + ' Variance')
    plt.ylabel('Distribution')
    plt.xlim(0, xmax)
    plt.title(gene+'_'+label)
    return fig  # save_figure(fig, gene+'_'+label, 'figures/permutationTests')


def plot_eg_hists(df, adata_obs, num_shuffles, gene, label):
    """ Plots the histogram of bootstrapped mean variances amongst clones vs. the true mean variance amongst clones"""
    """ Could another way to do this test be a KS test on all of the null distributions of mean variances generated vs.
     the true distribution of mean variances?"""
    labels_testing = adata_obs.copy()
    tested_gene = []
    statistics = []
    # Select the gene to test
    labels_testing[gene] = df.loc[:, gene]
    # get the mean variance of the true clonal groupings
    observedlabel_var = labels_testing.groupby(label)[gene].var()
    mean_observedlabel_variance = observedlabel_var.mean()
    # Generate a null distribution of mean variances by shuffling clonal groupings
    mean_shuffled_variances = []
    for i in range(num_shuffles):
        # create copy, honestly don't think this is necessary, but only hurts performance
        labels_testing_copy = labels_testing.copy(deep=True)
        # shuffle labels
        labels_testing_copy.loc[:,label] = np.random.permutation(labels_testing_copy[label].values)
        # Calculate variances of shuffled clonal groups
        shuffled_variances = labels_testing_copy.groupby(label)[gene].var()
        # Calculate mean of variances of shuffled clonal groups
        # mean_shuffled_variances.append(shuffled_variances.mean())
    # Plot
    fig, ax = plt.subplots(1, 1)
    data_shuffled = pd.Series(shuffled_variances)
    data = observedlabel_var
    xmax = data_shuffled.max() + 0.2
    bins = np.linspace(0, xmax, 20)
    plt.hist(data_shuffled, bins=bins, color='peru', alpha=0.5)
    plt.hist(data_shuffled, bins=bins, color='peru', histtype='step')
    plt.hist(data, bins=bins, color='midnightblue', alpha=0.5)
    plt.hist(data, bins=bins, color='midnightblue', histtype='step')
    # Plot Means for reference
    plt.axvline(data_shuffled.mean(), 0, 1, c='peru', ls='--')
    plt.axvline(data.mean(), 0, 1, c='midnightblue', ls='--')
    plt.yscale('log')
    plt.xscale('linear')
    plt.xlim(0, xmax)
    plt.xlabel('Variance within group')
    plt.ylabel('Events')
    plt.title(gene+'_'+label)
    plt.legend(['shuffled', 'true labels'])
    return fig, data, data_shuffled