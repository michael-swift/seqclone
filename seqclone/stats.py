# config parsing (.ini file)
import configparser
# standard data manipulation
import pandas as pd
import numpy as np
import sys
import os
# working with h5ad file
import scanpy as sc
# null distribution fitting 
from scipy.stats import norm
# bonferroni correction
from statsmodels.stats.multitest import multipletests

np.seterr(all='warn')

cfgFile = sys.argv[1]


def read_config(cfgFile):
    config = configparser.ConfigParser()
    config.read(cfgFile)
    stat_parameters = config['stat_parameters']
    io = config['IO']
    return stat_parameters, io, config


# Load the data and filter into a usable form
def prepare_data(CountsFile, datatype, highly_variable, n_highly_variable, onlyClones, remove_immune_receptors,
                normalize, filterCells):
    """ Accepts: H5ad file from scanpy where the adata.obs has a column "CLONE" denoting the clonal membership of the cell
        dataype: "scaled" or anything else would make it return log
        the rest of the parameters a passed by the config
    Returns: adata after filtering"""
    adata = sc.read_h5ad(CountsFile)
    # Use Scanpy to apply filtering criteria specified in config file
    adata, df = preprocess_wscanpy(adata, datatype, highly_variable, n_highly_variable, remove_immune_receptors,
                                  normalize, filterCells)
    # After filtering, select only cells which are clones
    if onlyClones:
        # Logic for dropping non-clones from the klein dataset
        # adata.obs.CLONE.fillna('None', inplace = True)
        adata = adata[adata.obs.CLONE != 'NaN', :]
        #Select only clones (applies to my dataset mostly)
        selector = adata.obs.CLONE.value_counts() > 1
        selector = selector[selector == True]
        adata = adata[adata.obs.CLONE.isin(selector.index), :]
        df = df[df.index.isin(adata.obs.index)]
   
    return adata, df


def preprocess_wscanpy(adata, datatype, highly_variable, n_highly_variable, remove_immune_receptors, normalize,
                      filterCells):
    if remove_immune_receptors:
        # TODO: somehow incorporate this information into the adata.var e.g immune-variable_receptor (True/False)
        immune_receptors = pd.read_csv('/metadata/immune_receptor_genes_keepConstantRegion.csv', index_col=0)
        immune_receptors.columns = ['genes']
        print("removing variable immune receptor genes which may drive clustering")
        adata = adata[:, ~adata.var.index.isin(immune_receptors.genes)]
    if filterCells:
        # Quality Filter Cells and Genes
        sc.pp.filter_cells(adata, min_genes=600, inplace=True)
        sc.pp.filter_cells(adata, min_counts=100000, inplace=True)
    # Filter out the lowest expressed genes for computation time
    sc.pp.filter_genes(adata, min_cells=50, inplace=True)
    sc.pp.filter_genes(adata, min_counts=200, inplace=True)
    print(adata.obs.shape, adata.var.shape, "shape of adata after filtering ")
    # Make parameter in cfg
    if normalize:
        sc.pp.normalize_total(adata, target_sum=1e6)
    sc.pp.log1p(adata, base=10)
    adata.raw = adata
    sc.pp.highly_variable_genes(adata, n_top_genes=n_highly_variable)
    # datatype logic
    if datatype == 'scaled':
        sc.pp.scale(adata)
    else:
        pass
    # Subset to highly variable gene
    if highly_variable:
        adata = adata[:, adata.var['highly_variable'] == True]
        highly_variable_genes = adata.var.index[adata.var["highly_variable"] == True]
    df = convert_sparse_to_dataframe(adata)
    return adata, df


def convert_sparse_to_dataframe(adata):
    """ Input: anndata object with sparse matrix as .X attribute
        Returns: Pandas dataframe with rows as cells and columns as genes
        My take: Converting the data to a df inefficient but convenient and possibly more readable"""
    # Get the gene expression values for each cell x gene
    columns = adata.var.index.to_list()
    index = adata.obs.index.to_list()
    try:
        dense_array = adata.X.toarray()
    except:
        dense_array = adata.X
    df = pd.DataFrame(data=dense_array, index=index, columns=columns)
    return df


def compare_variances(df, _adata_obs, num_shuffles, label, gene, well):
    """ comparer the mean variances of a shuffled labeling to the observed labeling of a given labeling and gene"""
    _adata_obs.loc[:, 'gene_name'] = df[gene]
    mean_shuffled_variances = []
    observedlabel_var = _adata_obs.groupby(label)['gene_name'].var()
    mean_observedlabel_variance = observedlabel_var.mean()
    
    
    # do the shuffling
    for i in range(num_shuffles):
        # create copy
        _adata_obs_copy = _adata_obs.copy(deep=True)
        # shuffle labels
        well = 'Experimental_Label'
        label = 'CLONE'
        list_of_dfs = []
        for group, frame in _adata_obs_copy.groupby(well):
            frame.loc[:,label] = np.random.permutation(frame.loc[:,label].values)
            list_of_dfs.append(frame)
        _adata_obs_copy = pd.concat(list_of_dfs)
        # groupby by label and compute variance
        shuffled_variances = _adata_obs_copy.groupby(label)['gene_name'].var()
        # Mean variance of every labeled group
        mean_shuffled_variances.append(shuffled_variances.mean())
    # make list into series TODO refactor to just add to a series?
    # This is the distribution of shuffled variances
    mean_shuffled_variances = pd.Series(mean_shuffled_variances)
    # Number of times shuffled variances are less than observed label variance,
    # higher number would be intragroup variance is higher
    test = mean_shuffled_variances <= mean_observedlabel_variance
    # less equal to observed (i.e. True's) by the number of tests
    # stat of 1 would be that shuffled variances always less or equal, 0 would be shuffled variances always more
    # this is a frequentist p value? kinda ... we'll call it a score
    gene_score = test.sum() / test.shape[0]
    return gene_score, gene, mean_shuffled_variances, mean_observedlabel_variance


def calculate_pvalue(mean_shuffled_variances, mean_observedlabel_variance):
    # Fit a normal distribution to the null variances I calculated
    mu, sigma = norm.fit(mean_shuffled_variances)
    if sigma == 0:
        # print('no p value possible for because null is not gaussian (sigma of zero)')
        pvalue = np.nan
    else:
        pvalue = norm.cdf(mean_observedlabel_variance, loc=mu, scale=sigma)
    return pvalue


def permutation_test(df, _adata_obs, num_shuffles, label, well):
    """ df is the cell x gene dataframe, adata_obs is metadataframe that
    contains cells as the index and a categorical column with the same name as the label"""
    # Get Annotation (clone data) and only what is in the scaled or transformed gene expression df
    _adata_obs = _adata_obs[_adata_obs.index.isin(df.index)]
    tested_genes = []
    gene_scores = []
    pvals = []
    genes = df.columns
    print("Running Permutation Test with", label, "and", num_shuffles, 'Shuffles')
    # This would be the natural place to parallelize:
    # Start an individual process which run the code in this loop a chunk of the genes iterable
    # merge all the results of the processes together after
    for gene in genes:
        gene_score, gene, mean_shuffled_variances, mean_observedlabel_variance = compare_variances(df, _adata_obs,
                                                                                                   num_shuffles, label,
                                                                                                   gene, well)
        tested_genes.append(gene)
        gene_scores.append(gene_score)
        # Calculate p-value by fitting Gaussian to null distribution
        # model = 'Gaussian' TODO could be other models
        pval = calculate_pvalue(mean_shuffled_variances, mean_observedlabel_variance)
        pvals.append(pval)
        print('Testing', gene)
    
    # merge results here
    scores_column = pd.Series(gene_scores)
    genes_column = pd.Series(tested_genes)
    pval_column = pd.Series(pvals)
    bonferroni_corrected_pvals = pd.Series(pvals)  # dummy column! not corrected yet!
    df_tests = pd.DataFrame([scores_column, genes_column, pval_column, bonferroni_corrected_pvals]).T
    df_tests.columns = ['score', 'gene', 'pvalue', 'corrected_pvalue']
    # what are the nas? probably where there is only 1 clone?
    df_tests.dropna(inplace=True)
    # multiple testing correction
    df_tests.loc[:, 'corrected_pvalue'] = multipletests(df_tests['pvalue'].values, alpha=0.05, method='bonferroni')[1]
    return df_tests


def write_results(df_tests, parameters, io, adata, config):
    # Merge with adata.var to get gene qc information
    df_results = pd.merge(adata.var, df_tests, left_index=True, right_on='gene')
    df_results.sort_values('corrected_pvalue', inplace=True)
    out_dir = io['out_dir']
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    df_results.to_csv('{}/{}_{}.csv'.format(io['out_dir'], parameters['label'], 'tests'))
    with open('{}/{}_configs.ini'.format(io['out_dir'], parameters['label']), 'w') as configfile:
        config.write(configfile)
    return df_results


# Main Function
def main():
    print("Running Stats on Variance Within Clones!")
    # Read the config file for parameters and input output information
    parameters, io, config = read_config(cfgFile)
    print(io['CountsFile'], "Using this data")
    # Apply filtering criteria for the data
    adata, df = prepare_data(io['CountsFile'],
                             parameters['datatype'],
                             parameters.getboolean('highly_variable'),
                             int(parameters['n_highly_variable']),
                             parameters.getboolean('onlyClones'),
                             parameters.getboolean('remove_immune_receptors'),
                             parameters.getboolean('normalize'),
                             parameters.getboolean('filterCells'))
    # Run the Tests
    print(adata)
    df_tests = permutation_test(df, adata.obs, int(parameters['num_shuffles']), parameters['label'], parameters['well'])
    # Write test results out as well as the exact config file used for this run (for logging)
    write_results(df_tests, parameters, io, adata, config)
    print("Completed Stats Testing... \n Good Job Everyone")
    return


if __name__ == "__main__":
    main()
