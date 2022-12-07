# config parsing (.ini file)
import configparser
# argparsing
import argparse
# standard data manipulation
import pandas as pd
import numpy as np
import os
# working with h5ad file
import scanpy as sc
# working with h5mu file
# null distribution fitting 
from scipy.stats import norm
np.seterr(all='warn')
parser = argparse.ArgumentParser(description='Runs permutation tests on h5ad formatted Adata object')
parser.add_argument('-config', type=str, default = 'config',
                    help='path to config file containing parameters for filtering')
parser.add_argument('-outdir', type = str, default = '.', 
                    help = 'output directory '
                            '(default: working directory)')
args = parser.parse_args()

## parse i/o paths
cfgFile = args.config
outdir = args.outdir

def read_config(cfgFile):
    config = configparser.ConfigParser()
    config.read(cfgFile)
    stat_parameters = config['stat_parameters']
    io = config['IO']
    return stat_parameters, io, config

# Load the data and filter into a usable form
def prepare_data(anndata, datatype, onlyclones, log,
                normalize, filterCells, layer):
    """ Accepts: H5ad file from scanpy where the adata.obs has a column "clone_id" denoting the clonal membership of the cell
        dataype: "scaled" or anything else would make it return log
        the rest of the parameters a passed by the config
    Returns: adata after filtering"""
    adata = sc.read_h5ad(anndata)
    # Use Scanpy to apply filtering criteria specified in config file
    adata, df = preprocess_wscanpy(adata, datatype, onlyclones,
                                   normalize, log, filterCells, layer)
    # After filtering, select only cells which are clones
    if onlyclones:
        adata = adata[adata.obs.clone_id != 'NaN', :]
        #Select only clones (applies to my dataset mostly)
        selector = adata.obs.clone_id.value_counts() > 1
        selector = selector[selector == True]
        adata = adata[adata.obs.clone_id.isin(selector.index), :]
        df = df[df.index.isin(adata.obs.index)]
    # Calculate qc metrics only after all filtering is done
    sc.pp.calculate_qc_metrics(adata, inplace=True)
    print("{} cells after filtering and {} genes after filtering".format(adata.obs.shape[0], adata.var.shape[0]))
    return adata, df


def preprocess_wscanpy(adata, datatype, onlyclones,
                       normalize, log, filterCells, layer):
    adata.raw = None
    if layer == 'raw_counts':
        adata.X = adata.layers['raw_counts']
    if filterCells:
        # Quality Filter Cells and Genes
        sc.pp.filter_cells(adata, min_genes=400, inplace=True)
        sc.pp.filter_cells(adata, min_counts=200, inplace=True)
    # Filter out the lowest expressed genes for computation time
    sc.pp.filter_genes(adata, min_cells=10, inplace=True)
    sc.pp.filter_genes(adata, min_counts=10, inplace=True)
    if normalize:
        sc.pp.normalize_total(adata, target_sum=1e4)
    if log:
        sc.pp.log1p(adata, base=2)
    adata.raw = adata
    #sc.pp.highly_variable_genes(adata)
    # datatype logic
    if datatype == 'scaled':
        sc.pp.scale(adata)
    df = convert_sparse_to_dataframe(adata)
    return adata, df

def shuffleWithinGroup(adata, label, group_name):
    """ returns a shuffled series of the label, shuffled within the groups supplied """
    list_of_dfs = []
    for group, frame in adata.obs.groupby(group_name):
        mframe.loc[:,label] = np.random.permutation(frame.loc[:,label])
        list_of_dfs.append(frame)
    df = pd.concat(list_of_dfs)
    return df.loc[:,label]

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

def calculateStatsDf(adata, label, well, stat, shuffle_first):
    df = convert_sparse_to_dataframe(adata)
    # Create datastructure to put test results
    stats = pd.DataFrame(df.columns)
    stats.set_index(0, inplace=True)
    # compute true statistic first
    if shuffle_first == True:
        adata.obs.loc[:,label] = shuffleWithinGroup(adata, label, well)
    else:
        df.loc[:,label] = adata.obs.clone_id
    df.loc[:,label] = adata.obs.clone_id
    # calculate variance of every get within shuffled clonal groups
    if stat == 'std':
        df_stats_shuffled = df.groupby(label).std()
    if stat == 'var':
        df_stats_shuffled = df.groupby(label).var()
    # calculate the mean variance of each gene when this permutation of labels is used
    means = df_stats_shuffled.mean().rename('true')
    stats = pd.merge(stats, means, left_index=True, right_index=True)
    return stats

def calculateShuffledParameter(adata, df, stats, num_shuffles, label, well, stat):
    modulus = num_shuffles / 100
    for i in range(num_shuffles):
        if i % modulus == 0:
            print("{}th shuffle".format(i))
        # perform a shuffle, only shuffling within the defined group e.g. lane / donor / timepoint
        adata.obs['shuffled'] = shuffleWithinGroup(adata, label, well)
        # add new shuffled column to df
        df.loc[:,label] = adata.obs['shuffled']
        # calculate variance of every get within shuffled clonal groups
        if stat == 'std':
            df_stats_shuffled = df.groupby(label).std()
        if stat == 'var':
            df_stats_shuffled = df.groupby(label).var()
        # calculate the mean variance of each gene when this permutation of labels is used
        means = df_stats_shuffled.mean().rename('null_'+ str(i))
        stats = pd.merge(stats, means, left_index=True, right_index=True)
    return stats

def _calculate_pvalue(mean_shuffled_variances, mean_observedlabel_variance):
    # Fit a normal distribution to the null variances I calculated
    mu, sigma = norm.fit(mean_shuffled_variances)
    if sigma == 0:
        # no p value possible for because null is not gaussian (sigma of zero)'
        pvalue = np.nan
    else:
        pvalue = norm.cdf(mean_observedlabel_variance, loc=mu, scale=sigma)
    return pvalue

def _createSummarydf(stats, adata):
    print("calculating p-values and effect sizes")
    frequentist_ps = []
    parametric_p = []
    effect_size = []
    for index, row in stats.iterrows():
        tests = row[1:]
        true = row[0]
        count = 0
        for i in tests:
            if i < true:
                count +=1
        frequentist_ps.append(count)
         
        parametric_p.append(_calculate_pvalue(tests, true))
        effect_size.append(true - np.mean(tests))
    summarydf = pd.DataFrame(frequentist_ps)
    summarydf.index = stats.index
    summarydf.columns = ['freq_pvalue']
    summarydf.loc[:,'pvalue'] = parametric_p
    summarydf.loc[:,'effect_size'] = effect_size
    summarydf = pd.merge(summarydf, adata.var, left_index=True, right_index=True)
    summarydf['clonally_expressed'] = summarydf['pvalue'] < 0.01
    summarydf.loc[:,'normalized_effect_size'] = summarydf['effect_size'] / summarydf['mean_counts']
    return summarydf


def _write_results(summarydf, stats, parameters, io, adata, config):
    # Merge with adata.var to get gene qc information
    summarydf.sort_values('pvalue', inplace=True)
    out_dir = io['out_dir']
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    summarydf.to_csv('{}/{}_{}.csv'.format(io['out_dir'], parameters['label'], 'var_df'))
    stats.to_csv('{}/{}_{}.csv'.format(io['out_dir'], parameters['label'], 'stats'))
    with open('{}/{}_configs.ini'.format(io['out_dir'], parameters['label']), 'w') as configfile:
        config.write(configfile)
    return summarydf


# Main Function
def main():
    print("Running Stats on Variance Within Clones!")
    # Read the config file for parameters and input output information
    parameters, io, config = read_config(cfgFile)
    print("Using this file: ", io['anndata'])
    # Prepare the anndata file
    adata, df = prepare_data(io['anndata'],
                             parameters['datatype'],
                             parameters.getboolean('onlyclones'),
                             parameters.getboolean('normalize'),
                             parameters.getboolean('log'),
                             parameters.getboolean('filterCells'), parameters['layer'])
    # Run the Tests
    # calculate True Variances for all genes within clones
    stats = calculateStatsDf(adata, parameters['label'], parameters['well'], parameters['stat'], parameters.getboolean('shuffle_first'))
    # calculate variances for permuted labels for all genes withing clones
    stats = calculateShuffledParameter(adata, df, stats, 
                                       int(parameters['num_shuffles']), 
                                       parameters['label'], 
                                       parameters['well'], parameters['stat'])
    summarydf = _createSummarydf(stats, adata)
    # Write test results out as well as the exact config file used for this run (for logging)
    _write_results(summarydf, stats, parameters, io, adata, config)
    print("Completed Stats Testing... \n Good Job Everyone")
    return


if __name__ == "__main__":
    main()
