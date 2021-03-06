# Import Libraries, some are uncessary right now

#config parsing (.ini file)
import configparser
#standard data manipulation
import pandas as pd
import numpy as np
import sys
import os
import random
import copy
import math
# working with h5ad file
import scanpy as sc
# plotting
from matplotlib import pyplot as plt
import matplotlib as mpl
import seaborn as sns
# null distribution fitting
from scipy.stats import norm
from scipy.stats import t
# bonferroni correction
from statsmodels.stats.multitest import multipletests

np.seterr(all = 'warn')

cfgFile = sys.argv[1] # '../switchy/SS2.ini'

def readConfig(cfgFile):
    config = configparser.ConfigParser()
    config.read(cfgFile)
    stat_parameters = config['stat_parameters']
    io = config['IO']
    CountsFile = io['CountsFile']
    out_dir = io['out_dir']
    return stat_parameters, io, config

# Load the data a get filter into a usable form
def prepareData(CountsFile, datatype, highly_variable, n_highly_variable, onlyClones, remove_immune_receptors, normalize, filterCells):
    """ Accepts: H5ad file from scanpy where the adata.obs has a column "CLONE" denoting the clonal membership of the cell
        dataype: "scaled" or anything else would make it return log
        the rest of the parameters a passed by the config
    Returns: adata after filtering"""
    adata = sc.read_h5ad(CountsFile)
    adata, df = preprocessWScanpy(adata, datatype, highly_variable, n_highly_variable, remove_immune_receptors, normalize, filterCells)
    # After filtering select only cells which are clones
    if onlyClones == True:
        # Logic for dropping non-clones from the klein dataset
        adata.obs.CLONE.fillna('None', inplace = True)
        adata = adata[adata.obs.CLONE != 'NaN' ,:]
        # Select only clones (applies to my dataset mostly)
        selector = adata.obs.CLONE.value_counts() > 1
        selector = selector[selector == True]
        adata = adata[adata.obs.CLONE.isin(selector.index), :]
        df = df[df.index.isin(adata.obs.index)]
    return adata, df


def preprocessWScanpy(adata, datatype, highly_variable, n_highly_variable, remove_immune_receptors, normalize, filterCells):
    if remove_immune_receptors == True:
        # TODO: Supply this as a metadata file or somehow incorporate this information into the adata.var e.g immune-variable_receptor (True/False)
        immune_receptors = pd.read_csv('/home/mswift/B_cells/CSR/sc_RNAseq/data_tables/metadata/immune_receptor_genes_keepConstantRegion.csv', index_col=0)
        immune_receptors.columns = ['genes']
        print("removing variable immune receptor genes which may drive clustering")
        adata = adata[:, ~adata.var.index.isin(immune_receptors.genes)]
    if filterCells == True:
    # Quality Filter Cells and Genes
        sc.pp.filter_cells(adata, min_genes=600, inplace = True)
        sc.pp.filter_cells(adata, min_counts=100000, inplace = True)
    # Filter out the lowest expressed genes for computation time
    sc.pp.filter_genes(adata, min_cells=50, inplace = True)
    sc.pp.filter_genes(adata, min_counts=200, inplace = True)
    print(adata.obs.shape, adata.var.shape, "shape of adata after filtering ")
    # Make parameter in cfg
    if normalize == True:
        sc.pp.normalize_total(adata)
    sc.pp.log1p(adata, base = 10)
    adata.raw = adata
    sc.pp.highly_variable_genes(adata, n_top_genes=n_highly_variable)
    # datatype logic
    if datatype == 'scaled':
        sc.pp.scale(adata)
    else:
        pass
    #Subset to highly variable gene
    if highly_variable == True:
        adata = adata[:,adata.var['highly_variable'] == True]
        highly_variable_genes = adata.var.index[adata.var["highly_variable"] == True]
    df = convertSparsetoDataFrame(adata)
    return adata, df

def convertSparsetoDataFrame(adata):
    """ Input: anndata object with sparse matrix as .X attribute
        Returns: Pandas dataframe with rows as cells and columns as genes
        My take: Converting the data to a df inefficient but convenient and possibly more readable"""
    # Get the gene expression values for each cell x gene
    columns = adata.var.index.to_list()
    index = adata.obs.index.to_list()
    try:
        denseArray = adata.X.toarray()
    except:
        denseArray = adata.X
    df = pd.DataFrame(data = denseArray, index = index , columns = columns )
    return df


################
### Plotting ###
################

def plotWaterfall(df, adata_obs, gene, label, filterzeroes):
    """ Strip plot where x axis is rank ordered by mean gene expression within a clone"""
    LabelsTesting = adata_obs
    LabelsTesting.loc[:,gene] = df[gene]
    if filterzeroes == True:
        LabelsTesting = LabelsTesting[LabelsTesting[gene] > 0]
        print(LabelsTesting.shape)
        selector = LabelsTesting[label].value_counts() > 1
        selector = selector[selector == True].index
        LabelsTesting = LabelsTesting[LabelsTesting[label].isin(selector)]
        print(LabelsTesting.shape)
    fig, ax = plt.subplots(1,1)
    order = LabelsTesting.groupby(label)[gene].mean().sort_values(ascending = False).index
    g = sns.stripplot(ax=ax, data = LabelsTesting, x = LabelsTesting[label], y = gene, order = order, color = None)
    #save_figure(fig, '{}_{}'.format(label, str(gene)))
    return g, LabelsTesting

def plotSwarm(df, adata_obs, gene, label):
    """ Strip plot where x axis is rank ordered by mean gene expression within a clone"""
    LabelsTesting = adata_obs
    fig, ax = plt.subplots(1,1)
    LabelsTesting.loc[:,gene] = df[gene]
    order = LabelsTesting.groupby(label)[gene].mean().sort_values(ascending = False).index
    g = sns.swarmplot(ax=ax, data = LabelsTesting, x = label, y = gene, order = order, color = None)
    #save_figure(fig, '{}_{}'.format(label, str(gene)))
    return g

def plotCI(df, adata, num_shuffles, gene, label, alpha):
    # This is time consuming to do twice, like this because I really am only plotting a subset of hits
    LabelsTesting = pd.merge(adata.obs[label], df[gene], left_index=True, right_index=True)
    tested_gene = []
    statistics = []
    # get the ordered means of the true labeling
    observedlabel_mean = LabelsTesting.groupby(label)[gene].mean().sort_values(ascending = False)
    ci_df = pd.DataFrame(observedlabel_mean)
    # shuffling loop
    for i in range(num_shuffles):
        # create copy out of superstition
        LabelsTestingCopy = LabelsTesting.copy(deep = True)
        # shuffle labels
        LabelsTestingCopy[label] = np.random.permutation(LabelsTestingCopy[label].values)
        # calculate mean gene expression of the shuffled
        shuffled_means = LabelsTestingCopy.groupby(label)[gene].mean()
        # add the shuffled mean to the ci_df
        ci_df = pd.merge(ci_df, shuffled_means, left_index=True, right_index=True)
    # First column is the true label means
    true_ordered_means = ci_df.iloc[:,0]
    # rest of the columns are the shuffled means
    shuffled_means = ci_df.iloc[:,1:]
    # Using T distribution
    T_low = shuffled_means.apply(lambda row: t.interval(alpha, row.shape[0]-1, loc = row.median(), scale=row.sem())[0], axis = 1)
    T_high = shuffled_means.apply(lambda row: t.interval(alpha, row.shape[0]-1, loc = row.median(), scale=row.sem())[1], axis = 1)
    lower_quantile = shuffled_means.quantile(q = 0.025, axis = 1)
    upper_quantile = shuffled_means.quantile(q = 0.975, axis = 1)

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

def plotTestHist(df, adata_obs, num_shuffles, gene, label, filterzeroes):
    """ Plots the histogram of bootstrapped mean variances amongst clones vs. the true mean variance amongst clones""" 
    """ Could another way to do this test be a KS test on all of the null distributions of mean variances generated vs. the true distribution of mean variances"""
    LabelsTesting = adata_obs.copy()
    tested_gene = []
    statistics = []
    # Select the gene to test
    LabelsTesting[gene] = df.loc[:,gene]
    if filterzeroes == True:
        LabelsTesting = LabelsTesting[LabelsTesting[gene] > 0]
        print(LabelsTesting.shape)
        selector = LabelsTesting[label].value_counts() > 1
        selector = selector[selector == True].index
        LabelsTesting = LabelsTesting[LabelsTesting[label].isin(selector)]
        print(LabelsTesting.shape)
    # get the mean variance of the true clonal groupings
    observedlabel_var = LabelsTesting.groupby(label)[gene].var()
    mean_observedlabel_variance = observedlabel_var.mean()
    # Generate a null distribution of mean variances by shuffling clonal groupings
    mean_shuffled_variances = []
    for i in range(num_shuffles):
        # create copy, honestly don't think this is necessary, but only hurts performance
        LabelsTestingCopy = LabelsTesting.copy(deep = True)
        # shuffle labels
        LabelsTestingCopy.loc[:,label] = np.random.permutation(LabelsTestingCopy[label].values)
        # Calculate variances of shuffled clonal groups
        shuffled_variances = LabelsTestingCopy.groupby(label)[gene].var()
        # Calculate mean of variances of shuffled clonal groups
        mean_shuffled_variances.append(shuffled_variances.mean())
    # Plot
    fig, ax = plt.subplots(1,1)
    data = pd.Series(mean_shuffled_variances)
    xmax = data.max() + 0.2
    bins = np.linspace(0, xmax, 100)
    plt.hist(data, bins = bins, color = 'peru', alpha = 0.5)
    plt.hist(data, bins = bins, color = 'peru', histtype='step')
    plt.axvline(mean_observedlabel_variance, 0, 1, c = 'midnightblue', ls = '--')
    plt.yscale('log')
    plt.xscale('linear')
    plt.xlabel('Average ' + label + ' Variance')
    plt.ylabel('Distribution')
    plt.xlim(0, xmax)
    plt.title(gene+'_'+label)
    return fig#save_figure(fig, gene+'_'+label, 'figures/permutationTests')

def plotEgHists(df, adata_obs, num_shuffles, gene, label):
    """ Plots the histogram of bootstrapped mean variances amongst clones vs. the true mean variance amongst clones""" 
    """ Could another way to do this test be a KS test on all of the null distributions of mean variances generated vs. the true distribution of mean variances?"""
    LabelsTesting = adata_obs.copy()
    tested_gene = []
    statistics = []
    # Select the gene to test
    LabelsTesting[gene] = df.loc[:,gene]
    # get the mean variance of the true clonal groupings
    observedlabel_var = LabelsTesting.groupby(label)[gene].var()
    mean_observedlabel_variance = observedlabel_var.mean()
    # Generate a null distribution of mean variances by shuffling clonal groupings
    mean_shuffled_variances = []
    for i in range(num_shuffles):
        # create copy, honestly don't think this is necessary, but only hurts performance
        LabelsTestingCopy = LabelsTesting.copy(deep = True)
        # shuffle labels
        LabelsTestingCopy.loc[:,label] = np.random.permutation(LabelsTestingCopy[label].values)
        # Calculate variances of shuffled clonal groups
        shuffled_variances = LabelsTestingCopy.groupby(label)[gene].var()
        # Calculate mean of variances of shuffled clonal groups
        #mean_shuffled_variances.append(shuffled_variances.mean())
    # Plot 
    fig, ax = plt.subplots(1,1)
    data_shuffled = pd.Series(shuffled_variances)
    data = observedlabel_var
    xmax = data_shuffled.max() + 0.2
    bins = np.linspace(0, xmax, 20)
    plt.hist(data_shuffled, bins = bins, color = 'peru', alpha = 0.5)
    plt.hist(data_shuffled, bins = bins, color = 'peru', histtype='step')
    plt.hist(data, bins = bins, color = 'midnightblue', alpha = 0.5)
    plt.hist(data, bins = bins, color = 'midnightblue', histtype='step')
    # Plot Means for reference
    plt.axvline(data_shuffled.mean(), 0, 1, c = 'peru', ls = '--')
    plt.axvline(data.mean(), 0, 1, c = 'midnightblue', ls = '--')
    plt.yscale('log')
    plt.xscale('linear')
    plt.xlim(0, xmax)
    plt.xlabel('Variance within group')
    plt.ylabel('Events')
    plt.title(gene+'_'+label)
    plt.legend(['shuffled', 'true labels'])
    return fig, data, data_shuffled
    #save_figure(fig, gene+'_'+label, 'figures/permutationTests')

def compareVariances(df, LabelsTesting, num_shuffles, label, gene):
    "For each gene compare the mean variances of a shuffled labeling to the observed labeling"
    LabelsTesting.loc[:,'gene_name'] = df[gene]
    mean_shuffled_variances = []
    observedlabel_var = LabelsTesting.groupby(label)['gene_name'].var()
    mean_observedlabel_variance = observedlabel_var.mean()
    # do the shuffling
    for i in range(num_shuffles):
        # create copy
        LabelsTestingCopy = LabelsTesting.copy(deep = True)
        #shuffle labels
        LabelsTestingCopy.loc[:, label] = np.random.permutation(LabelsTesting[label].values)
        # groupby by label and compute variance
        shuffled_variances = LabelsTestingCopy.groupby(label)['gene_name'].var()
        # Mean variance of every labeled group
        mean_shuffled_variances.append(shuffled_variances.mean())
    #make list into series TODO refactor to just add to a series?
    #This is the distribution of shuffled variances
    mean_shuffled_variances = pd.Series(mean_shuffled_variances)
    # Number of times shuffled variances are less than observed label variance, higher number would be intragroup variance is higher
    test = mean_shuffled_variances <= mean_observedlabel_variance
        # less equal to observed (i.e. True's) by the number of tests
        # stat of 1 would be that shuffled variances always less or equal, 0 would be shuffled variances always more 
        # this is a frequentist p value? kinda ... we'll call it a score
    gene_score = test.sum() / test.shape[0]
    return gene_score, gene, mean_shuffled_variances, mean_observedlabel_variance

def calculatePvalue(mean_shuffled_variances, mean_observedlabel_variance):
    # Fit a normal distribution to the null variances I calculated
    mu, sigma = norm.fit(mean_shuffled_variances)
    if sigma == 0:
        #print('no p value possible for because null is not gaussian (sigma of zero)')
        pvalue = np.nan
    else:
        pvalue = norm.cdf(mean_observedlabel_variance, loc = mu, scale = sigma)
    return pvalue

def permutationTest(df, adata_obs, num_shuffles, label):
    """ df is the cell x gene dataframe, adata_obs is metadataframe that
    contains cells as the index and a categorical column with the same name as the label"""
    # Get Annotation (clone data) and only what is in the scaled or transformed gene expression df
    LabelsTesting = adata_obs[adata_obs.index.isin(df.index)]
    tested_genes = []
    gene_scores = []
    pvals = []
    genes = df.columns
    print("Running Permutation Test with", label, "and", num_shuffles, 'Shuffles')
    genes_tested = 0
    for gene in genes:
        gene_score, gene, mean_shuffled_variances, mean_observedlabel_variance = compareVariances(df, LabelsTesting, num_shuffles, label, gene)
        tested_genes.append(gene)
        gene_scores.append(gene_score)
        # Calculate p-value by fitting Gaussian to null distribution
        # model = 'Gaussian' TODO could be other models
        pval = calculatePvalue(mean_shuffled_variances, mean_observedlabel_variance)
        pvals.append(pval)
        print('Testing', gene)
        genes_tested += 1
        print(genes_tested, "genes tested thusfar")
    scoresColumn = pd.Series(gene_scores)
    genesColumn = pd.Series(tested_genes)
    pvalColumn = pd.Series(pvals)
    bonferroniCorrected_pvals = pd.Series(pvals) # dummy column! not corrected yet!
    df_tests = pd.DataFrame([scoresColumn, genesColumn, pvalColumn, bonferroniCorrected_pvals]).T
    df_tests.columns = ['score', 'gene', 'pvalue', 'corrected_pvalue']
    # what are the nas? probably where there is only 1 clone?
    df_tests.dropna(inplace = True)
    try:
        df_tests.loc[:,'corrected_pvalue'] = multipletests(df_tests['pvalue'].values, alpha = 0.05, method = 'bonferroni')[1]
    except:
        df_tests.loc[:,'corrected_pvalue'] = "NotCorrected"
    return df_tests

def write_results(df_tests, parameters, io, adata, config):
    df_results = pd.merge(adata.var, df_tests, left_index = True, right_on='gene')
    df_results.sort_values('corrected_pvalue', inplace = True)
    out_dir = io['out_dir']
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    df_results.to_csv('{}/{}_{}.csv'.format(io['out_dir'], parameters['label'], 'tests'))
    with open('{}/{}_configs.ini'.format(io['out_dir'], parameters['label']), 'w') as configfile:
        config.write(configfile)
    return df_results


###################
## Main Function ## 
def main():
    print("Running Stats on Variance Within Clones!")
    # Read the config file for parameters and input output information
    parameters, io, config = readConfig(cfgFile)
    print(io['CountsFile'], "Using this data")
    # Apply filtering criteria for the data
    adata, df = prepareData(io['CountsFile'],
                                parameters['datatype'],
                                parameters.getboolean('highly_variable'),
                                int(parameters['n_highly_variable']),
                                parameters.getboolean('onlyClones'),
                                parameters.getboolean('remove_immune_receptors'),
                                parameters.getboolean('normalize'),
                                parameters.getboolean('filterCells'))
    # Run the Tests
    df_tests = permutationTest(df, adata.obs, int(parameters['num_shuffles']), parameters['label'])
    # Write test results out as well as the exact config file used for this run (for logging)
    df_tests = write_results(df_tests, parameters, io, adata, config)
    print("Completed Stats Testing... \n Good Job Everyone")
    return

if __name__ == "__main__":
    main()
