import pandas as pd
import scipy
import numpy as np
import seaborn as sns
import matplotlib as mpl
#from sinaplot import sinaplot
import scanpy as sc
from matplotlib import pyplot as plt
import scanpy.external as sce
import os 
import scipy.spatial.distance
cwd = os.getcwd()
print(cwd)

def getLineagesFromChangeo(changeodb, print_summary):
    """subsets the changeo_db output by bracer by only those cells which are within lineages (non singletons)"""
    df = changeodb
    _df = df[df.CLONE != "None"] # get rid of unassigned cells (no BCR reconstructed)
    _df = (df.CLONE.value_counts() > 1) #find clones with more than 1 member
    if print_summary == True:
        print( "There are", len(_df[_df == 1]), "lineages with more than one member")
    CHANGEO_confidentlineages = df[df.CLONE.isin(_df[_df == 1].index)].sort_values('CLONE')
    CHANGEO_confidentlineages = CHANGEO_confidentlineages[CHANGEO_confidentlineages.CLONE != 'None']
    if print_summary == True:
        print("number of cells in original dataframe", df.shape[0])
        print("number of distinct Clones in original dataframe", df.drop_duplicates('CLONE').shape[0] -1) #subtract 1 for the 'None' entry
        print(CHANGEO_confidentlineages.shape[0]/df.shape[0], 'percent of cells in a lineage' )
    return CHANGEO_confidentlineages

# calculate distance metric
def calculate_distance_metric(_clonal_comparison, df):
    # iterate, i know this is bad style.. 
    for index, row in _clonal_comparison.iterrows():
        # Cell to cell comparison
        cell1 = row['cell1']
        cell2 = row['cell2']
        # df.loc is an array that you make
        _clonal_comparison.loc[index, 'distance'] = distance.euclidean(df.loc[cell2, :].values, df.loc[cell1, :].values)
    return _clonal_comparison

# get lineage correlations for all groups
def calculateLineageCorrelations(counts_table, changeo_db, method):
    """ returns within lineage, amongst lineage, amongst all cells 
    ie a list or series of correlations between different groupings of cells (all, sisters, within lineage etc.)"""

    if method == 'euclidean_distance':
        method = scipy.spatial.distance.euclidean
    lineage, within_lineage = calculateCorrelationsWithinLineage(changeo_db, counts_table, method)
    #TODO seemlessly incorporate the print summary logic (maybe i don't really even need it at this point)
    amongst_lineage = calculateCorrelations(counts_table[counts_table.index.isin(getLineagesFromChangeo(changeo_db, print_summary=False).CELL)],
                               method)
    amongst_all_cells = calculateCorrelations(counts_table, method)
    return within_lineage, amongst_lineage, amongst_all_cells

# lineage correlation functions

def calculateCorrelations(df, method):
    # create empty list
    correlations_list = []
    #transform matrix to take pearson correllation
    _pearsonCorrAll = df.T.corr(method=method)

    for i in range(len(_pearsonCorrAll)):
        correlations_list += list(_pearsonCorrAll.iloc[(i+1):,i])
    correlations = pd.Series(correlations_list).dropna()
    return correlations

def calculateCorrelationsWithinLineage(changeodb, counts_table, method):
    """idf is the dataframe used to index should be the changeo df, or minimally df with cell name and its clone information in columns,
    counts_table is the cell by gene counts table from which the vectors of comparison come"""
    idf = changeodb
    lineage_ziplist = []
    correlations = []
    #iterate through each unique lineage ID
    for clone in idf.CLONE.unique():
        #get series of each cell in the lineage
        series_of_cells = idf.CELL[idf.CLONE == clone]
        #if the lineage is larger than 2 do the comparison
        if len(series_of_cells.unique()) < 2:
            pass
        else:
            #print(clone)
            lineage_counts_table = counts_table[counts_table.index.isin(series_of_cells)]
            correlation_matrix = lineage_counts_table.T.corr(method=method)
            ## subset array where no duplicate comparisons and no 1's
            mask = np.zeros_like(correlation_matrix, dtype=np.bool)
            mask[np.triu_indices_from(mask)] = True
            drop_dup = correlation_matrix.where(mask)
            drop_dup = drop_dup.values.flatten()
            drop_dup = drop_dup[drop_dup != 1]
            drop_dup = drop_dup[~np.isnan(drop_dup)] #not clear when and why I was getting nan's
            #it may be when I filter the gene counts data frame for good cells some of the clones don't make it
            correlations.extend(list(drop_dup))
            lineage_ziplist.append(clone) 
    return lineage_ziplist, correlations

def makeHumanReadableSJTable(df_IGH, exon_start_dict):
    """Adds columns to the SJout file which aid in interpreting the various splice junctions accepts and SJout formatted file and a dictionary which maps coordinates to starting exons
       Returns: SJout dataframes with extra columns, 1 df for the ab transcription (J spliced to downstream constant region) and 1 with just the switch transcripts defined as transcripts where the splice donor to exon1 is not a J gene """
    all_exon_starts = exon_start_dict
    _constant_regionspliceacceptors = df_IGH.copy()

    _J_genes_start_dict = {105865406: "IGHJ1", 105863197: "IGHJ6", 105864586: 'IGHJ3', 105865198: 'IGHJ2', 105864214: 'IGHJ4', 105863813: 'IGHJ5'}
    _constant_regionspliceacceptors['J_exon'] = _constant_regionspliceacceptors['end']
    _constant_regionspliceacceptors['J_exon'] = _constant_regionspliceacceptors['J_exon'].map(_J_genes_start_dict).fillna(_constant_regionspliceacceptors.end)
    # adding exon column
    _constant_regionspliceacceptors['exon_start'] = _constant_regionspliceacceptors['start'].replace(all_exon_starts)
    # filter by more than 2 unique mapping reads
    _constant_regionspliceacceptors = _constant_regionspliceacceptors[_constant_regionspliceacceptors.unique_mapping > 2]
    _constant_regionspliceacceptors = _constant_regionspliceacceptors.fillna('none')
    #filter by more than 12 overhang (gets rid of E mapping artefact which TODO I should investigate
    #cast the column to string
    _constant_regionspliceacceptors['exon_start'] = _constant_regionspliceacceptors["exon_start"].astype(str)
    # splices to the J genes
    _J_tx_df = _constant_regionspliceacceptors[_constant_regionspliceacceptors['J_exon'].isin(_J_genes_start_dict.values())]
    #take only the first exon of each constant region, this way we know we are counting switch transcripts
    #filter out the splicing counts coming from the J genes to the first exon
    switch_tx_df = _constant_regionspliceacceptors[(_constant_regionspliceacceptors.exon_start.str.contains('exon1'))
                                & (_constant_regionspliceacceptors.end < 105863197) # splice donor comes downstream of last J gene
                                & (~_constant_regionspliceacceptors.exon_start.str.contains('3prime'))] #get rid of odd but possibly interesting splicing where the acceptor is the 3prime? This may be related to backwards recombination or antisense transcription?
    #filter out where J is the acceptor of the splice (possibly interesting but not interested right now)
    _J_tx_df = _J_tx_df[~_J_tx_df.exon_start.str.contains("IGHJ")]
    _J_tx_df['exon_simple'] = _J_tx_df['exon_start'].str.split('_', expand = True)[0]
    return _J_tx_df, switch_tx_df



def addMetadataToAnnData(adata, changeodb_H, ab_tx):
    _adata = adata.copy()
    #_adata.obs = _adata.obs.reset_index()
    #split cell names which includes superfluous information about the Donor ID
    cell_names = _adata.obs.index.str[:-2]
    #set loom index to cell names
    _adata.obs.index = cell_names
    #Call Isotypes
    df_isotype_calls, x = callIsotypeByBracerSJout(changeodb_H, ab_tx, plot = False)
    #select most abundant isotype
    df_isotype_calls = df_isotype_calls.sort_values('unique_mapping',
                                                    ascending=False).drop_duplicates(subset='cell')
    #make a cell column for merging dataframes
    _adata.obs['cell'] = _adata.obs.index
    #perform the merge to get isotype information in scanpy
    _adata.obs = pd.merge(_adata.obs, df_isotype_calls, left_on='cell', right_on='cell', how='left')
    # set index to cell
    _adata.obs.set_index('cell', inplace=True)
    #drop "cell" name from the index in df
    _adata.obs.index.name = None
    # make new changeodb_H column called cell to facilitate merging
    changeodb_H['cell'] = changeodb_H.loc[:,'CELL']
    #Filter to only functional assemblies
    changeodb_H = changeodb_H[changeodb_H.FUNCTIONAL == True]
    # Drop any duplicates which would be possible doublets
    changeodb_H = changeodb_H.drop_duplicates(subset='cell')
    #perform the same merging dance as before
    _adata.obs['cell'] = _adata.obs.index
    _adata.obs = pd.merge(_adata.obs, changeodb_H, left_on= 'cell', right_on = 'cell', how = 'left')
    _adata.obs.set_index('cell', inplace = True)
    _adata.obs.index.name = None
    #_adata.obs = _adata.obs.fillna('no_assembly')
    return _adata

def loadSJoutIGH(filename, metadata):
    df_sjout = pd.read_feather(filename)
    #filter dataframe to just the IGH locus
    print("filtering SJout to just IGH locus")
    df_IGH = df_sjout[(df_sjout['end'] > 105550000)]
    #load metadata about constant region exon coordinates and I exons
    df_exoncoordinates = pd.read_csv(metadata, header = None, names = ['exon', 'coordinate'] )
    #Apply STAR vs. ENSMBL indexing correction
    df_exoncoordinates.loc[~df_exoncoordinates.exon.str.contains('exon1'), 'coordinate'] = df_exoncoordinates.coordinate+1
    
    df_exoncoordinates.loc[df_exoncoordinates.exon.str.contains('IGHG3_exon1'), 'coordinate'] = df_exoncoordinates.coordinate+1


    #load metadata about constant region exon coordinates and I exons
    #df_exoncoordinates = pd.read_csv(metadata, header = None, names = ['exon', 'coordinate'] )

    zipper = zip(df_exoncoordinates.coordinate, df_exoncoordinates.exon)
    all_exon_starts = dict(list(zipper))
    print("making SJTable human readable")
    ab_tx , switch_tx = makeHumanReadableSJTable(df_IGH, all_exon_starts)
    return ab_tx, switch_tx

def loadChangeoDbH(filepath):
    changeo_db = pd.read_csv(filepath, index_col = 0)
    changeo_db = changeo_db[changeo_db.LOCUS == 'H']
    changeo_db["CLONE"] = changeo_db['MERGE_CLONE'] #MERGE CLONE column was created in data combination process to make the column unique
    return changeo_db

def loadData(SJout, loom_data, gene_counts, changeo_db):
    ab_tx, switch_tx = loadSJoutIGH(SJout)
    ###
    # Assemblies (Changeo)
    changeo_db_H = loadChangeoDbH(changeo_db)

    ###
    # Gene counts
    ###
    print("loading anndata")
    adata = 'placeholder'
    loom_adata = sc.read_loom(loom_data)
    adata = sc.read_h5ad(gene_counts)
    return  ab_tx, switch_tx, adata, loom_adata, changeo_db_H

def preprocessGeneCounts(gene_counts, num_counted_reads):

    print("filtering cells with less than", str(num_counted_reads), "counted reads")
    gene_counts_filtered = gene_counts[gene_counts.sum(axis=1) > num_counted_reads] # cell have more than 30000 mapped reads
    print("normalizing to counts per million")
    _CPM = gene_counts_filtered.div(gene_counts_filtered.sum(axis = 1), axis=0) * 1e6 # counts per million
    print("log base 2 transforming")
    _logCPM = np.log2(_CPM + 1) #log base 2p
    return _logCPM

def preprocessScanpy(adata, num_counted_reads, num_genes, min_cells, n_neighbors, num_highly_variable):
    """ returns a new adata object that has been processed """
    _adata = adata.copy()
    _adata.var_names_make_unique()
    print("making var_names unique")
    print("filtering cells with less than", num_counted_reads, "counted reads")
    sc.pp.filter_cells(_adata,min_counts=num_counted_reads)
    print("filtering cells with less than", num_genes, "genes detected")
    sc.pp.filter_cells(_adata,min_genes=num_genes)
    print("filtering genes detected in less than", min_cells)
    sc.pp.filter_genes(_adata,min_cells=min_cells)
    sc.pp.calculate_qc_metrics(_adata, inplace = True)
    print("normalizing by total counts per cell")
    sc.pp.normalize_total(_adata, exclude_highly_expressed=True)
    print("log transforming data")
    sc.pp.log1p(_adata, base=10)
    _adata.raw = _adata
    # Remove ERCCs which could drive clustering based on which batch was used or whether or not they were spiked in
    ERCCs = _adata.var.index[_adata.var.index.str.contains("ERCC-")].to_list()
    _adata = _adata[:, ~_adata.var.index.isin(ERCCs)]
    print("removed ERCC sequences from genes to cluster on")
    # Remove Immune Receptor Genes which could drive clustering
    immune_receptors = pd.read_csv('/home/mswift/B_cells/CSR/sc_RNAseq/data_tables/metadata/immune_receptor_genes_keepConstantRegion.csv', index_col=0)
    immune_receptors.columns = ['genes']
    print("removing variable immune receptor genes which may drive clustering")
    _adata = _adata[:, ~_adata.var.index.isin(immune_receptors.genes)]
    print("calculating highly variable genes")
    sc.pp.highly_variable_genes(_adata, n_top_genes = num_highly_variable)
    sc.pp.scale(_adata)
    print("calculating PCA")
    sc.pp.pca(_adata)
    neighbors = n_neighbors
    print("creating neighbors graph with", n_neighbors)
    ## TODO some batch Correction?
    sc.pp.neighbors(_adata, n_neighbors=neighbors)
    print('umapping and tsne-ing')
    sc.tl.umap(_adata)
    sc.tl.tsne(_adata)
    return _adata

def plotPointPlotLocus(IGH_locus_df, cell_list, color):
    """makes a point plot of the IGH locus where reach observation is a cell and an observation consists
    of counts for each of the genes at the IgH locus"""
    sns.set(style = "whitegrid", context = 'paper')
    IGH_locus_df = IGH_locus_df[IGH_locus_df.index.isin(cell_list)]
    ## Point Plot individual clones
    #Data munging
    point_plot_df = IGH_locus_df

    point_plot_df = point_plot_df.reset_index()
    
    point_plot_df = pd.melt(point_plot_df, value_vars = point_plot_df.columns[1:], id_vars = 'index')
    point_plot_df.columns = ['cell', 'exon', 'log CPM']
        
    # plotting
    fig, ax = plt.subplots(1, 1, figsize=(6,3))
    sns.set_palette(color)
    sns.pointplot(data = point_plot_df,
                  x = 'exon', y ='log CPM', hue='cell', dodge =.2,
                  join = True, legend=None, alpha = 0.5)
    ax.set_ylabel("log$_2$ CPM")
    ax.set_xlabel("")
    ax.legend_.remove()
    sns.pointplot(data = point_plot_df,
                  x = 'exon', y ='log CPM',
                  join = True, linestyles = "--", color = 'k')
    sns.despine()
    ax.set_ylabel("log$_2$ CPM")
    ax.set_xlabel("")
    sns.set_style("whitegrid", {'axes.grid' : False})
    return fig, ax

## Scanpy helper functions:
def getCellByGeneMatrix(adata):
    """input: adata
        return: df cell by gene"""
    # get matrix
    _x = pd.DataFrame(adata.X)
    # add column names (genes)
    _x.columns = adata.var_names
    # add row names (cells)
    _x = _x.set_index(adata.obs_names)
    return _x

def getEmbeddingCoordinates(adata, embedding):
    #filtering for good cells
    embedding_key = 'X_' + embedding
    _df_embedding = sc.get.obs_df(adata, obsm_keys = [(embedding_key, 0), (embedding_key, 1)])
    return _df_embedding

def plotClonalCorrelations(gene_counts, changeodb_H, method):
    """ Accepts gene counts data frame with cells as rows
    and genes of interest as columns, Changeodb, and the method """
    within_lineage, amongst_lineage, amongst_all_cells = calculateLineageCorrelations(gene_counts, changeodb_H, method)
    hist = plt.grid()
    if method == 'euclidean_distance':
        bins = np.linspace(0,40,30)
    else:
        bins = np.linspace(-1.1,1.1, 30)

    plt.hist(amongst_lineage, bins=bins, color = 'blue', label='Unrelated Pairs', density=True, histtype='stepfilled', alpha=0.5)
    #plt.hist(amongst_lineage, bins=bins, color = 'blue', label='Unrelated Pairs', density=True, histtype='step')
    # Pearson correlation hist for all pairs
    # Pearson correlation hist for related pairs
    plt.hist(within_lineage, bins=bins, color = 'red', label='Related Pairs', density=True, histtype='stepfilled', alpha = 0.5)
    #plt.hist(within_lineage, bins=bins, color = 'red', label='Related Pairs', density=True, histtype='step')
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.xlabel(method)
    plt.ylabel('Density')
    print("Mean Pearson correlation amongst all cells in grey:",amongst_all_cells.mean())
    print("Mean Pearson correlation within lineages in red:", pd.Series(within_lineage).mean(), 'consisting of', pd.Series(within_lineage).shape[0], "comparisons")
    print("Mean Pearson correlation amongst all lineages in blue:", (pd.Series(amongst_lineage).mean()), 'consisting of', pd.Series(amongst_lineage).shape[0], "comparisons")
    print(scipy.stats.ks_2samp(within_lineage, amongst_lineage), 'KS-test result comparing amongst all lineages to within lineages')


def plotLocusHeatmap(sj_out_df, isotype_call_df, row_colorby):
    """plots a clustermap using seaborn an sjout file and the isotype calls, rows can be colored by "condition" or by "isotype" """
    sj_out_df = sj_out_df[sj_out_df.cell.isin(isotype_call_df.cell)]
    #collects all switch transcripts to each constant region exon even if they come from different I-exon coordinates
    sj_out_df.loc[:,'cell'] = sj_out_df.loc[:,'cell'].str.replace('-','_')
    sj_out_df = sj_out_df[sj_out_df.exon_start.str.contains('exon1')]
    sum_df = sj_out_df.groupby(['cell', 'exon_start']).sum()
    # log transform the uniquely mapping reads
    sum_df['unique_mapping_log2'] = np.log2(sum_df['unique_mapping'])
    #put data in long form
    df = sum_df.unique_mapping_log2.unstack().fillna(0)
    # make a condition column (for a color bar)
    _df = df.copy()
    _df['condition'], _df['cell'] = _df.index.str.split('_', 1).str

    # map color dictionary to condition column
    if row_colorby == 'condition':
        # hard coded condition color bar based on plate names from the cluster
        cultureConditionDict = {'P1': 'firebrick', 'D2': 'firebrick', 'MS': 'firebrick', 'CTY':'firebrick', 'IL6811':'firebrick', 'NaiveBcellsplate1':'midnightblue', 'NaiveBcellsplate2':'midnightblue', 'PASL':'green', 'AgSeqPooled': 'midnightblue'}
        row_colors = _df.condition.map(cultureConditionDict)
    #for isotype coloring but could be about something else if necessary
    #TODO: still some Nans in the row_colors that I'm unsure of where they are coming from
    else:
        _df["CELL"] = _df['condition'] + "_" + _df['cell']
        isotypeColorDict = dict(zip(isotype_call_df['cell'], isotype_call_df['node_color']))
        row_colors = _df.CELL.map(isotypeColorDict)
    #plot correlogram
    # plot
    return sns.clustermap(df, row_colors=row_colors), _df


def callIsotypeBySJout(ab_tx, plot):
    #slow for large datasets
    # other approaches could be better 
    # e.g https://stackoverflow.com/questions/15705630/get-the-rows-which-have-the-max-value-in-groups-using-groupby
    _df = ab_tx.copy(deep=True)
    # Parse Exon Start column
    # find the J to 1st constant region junction with max support (max unique mapping reads), this could be decreasing sensitivity by not calling cells in the act
    idx = _df.groupby(['cell'], sort=False)['unique_mapping'].transform(max) == _df['unique_mapping']
    #filter sjout df by this
    isotype_calls_df = _df[idx]
    if plot == True:
        f, ax = plt.subplots(figsize=(6, 15))

        ax = sns.barplot(data=isotype_calls_df.exon_simple.value_counts().to_frame('counts').reset_index(), y = 'index', x = 'counts')
        ax.set(xlabel="Isotypes called by splice junctions")
        ax.set()
    isotype_calls_df.loc[:,'ISOTYPE_by_splice'] = isotype_calls_df.loc[:,'exon_start'].str.split('_', expand=True)[0]
    isotype_list = ['None', 'IGHM', 'IGHD', 'IGHG3', 'IGHDM', 'IGHG1', 'IGHA2', 'IGHG2', 'IGHG4', 'IGHE', 'IGHA1', 'nan']
    color_list = ['grey', 'green', 'green', 'red', 'green', 'black', 'blue', 'magenta', 'pink','cyan', 'blue', 'grey']
    color_isotype_dict = dict(zip(isotype_list, color_list))
    isotype_calls_df.loc[:,'node_color']= isotype_calls_df.ISOTYPE_by_splice.map(color_isotype_dict)
    return isotype_calls_df

def callIsotypeByBracer(changeodb):
    changeodb_H = changeodb[changeodb.LOCUS == 'H']
    f, ax = plt.subplots(figsize=(6, 15))

    ax = sns.barplot(data=changeodb_H.ISOTYPE.value_counts().to_frame('counts').reset_index(), y = 'index', x ='counts')
    ax.set(xlabel="Isotypes called by BrACER assembly")
    ax.set()
    return changeodb_H

def callIsotypeByBracerSJout(changeodb, ab_tx, plot):
    dfBracer = changeodb[changeodb.LOCUS == 'H'].copy()
    #only take functional Bracer assemblies (i.e. no stop codons thus the principal igh antibody being made)
    dfBracer = dfBracer[dfBracer.FUNCTIONAL == True]
    sj_out_df = ab_tx.copy()
    #make "CELL" column to facilitate merge with "cell"
    sj_out_df["CELL"] = sj_out_df.loc[:,"cell"]

    merged_df = pd.merge(dfBracer, sj_out_df, how='inner', on='CELL')
    #Take the the splice Junction calls where the J gene matches the J gene from the assembly
    _df = merged_df[merged_df.J_CALL == merged_df.J_exon][['J_exon', 'cell', 'ISOTYPE', 'exon_start', 'unique_mapping']] 
    # find the J to 1st constant region junction with max support (max unique mapping reads)
    idx = _df.groupby(['cell'], sort=False)['unique_mapping'].transform(max) == _df['unique_mapping']
    #filter sjout df by this
    isotype_calls_df = _df[idx]
    #Plot Bar plot of relative abundances of each isotype
    if plot == True:
        f, ax = plt.subplots(figsize=(6, 15))

        ax = sns.barplot(data=isotype_calls_df.exon_start.value_counts().to_frame('counts').reset_index(), y = 'index', x = 'counts')

        ax.set(xlabel="Isotypes called by combined Bracer and SJ.out")
        ax.set()
    else:
        ax = None
    # Simplify the exon_start column for plotting/display purposes
    isotype_calls_df.loc[:,'ISOTYPE_by_splice'] = isotype_calls_df.loc[:,'exon_start'].str.split('_', expand=True)[0]
    isotype_list = ['None', 'IGHM', 'IGHD', 'IGHG3', 'IGHDM', 'IGHG1', 'IGHA2', 'IGHG2', 'IGHG4', 'IGHE', 'IGHA1', 'nan']
    color_list = ['grey', 'green', 'green', 'red', 'green', 'black', 'blue', 'magenta', 'pink','cyan', 'blue', 'grey']
    color_isotype_dict = dict(zip(isotype_list, color_list))
    isotype_calls_df.loc[:,'node_color']= isotype_calls_df.ISOTYPE_by_splice.map(color_isotype_dict)
    return isotype_calls_df, ax

def findActofSwitchingCells(ab_sjoutdf, threshold, changeodb_H, plot):
    """ab_sjoutdf is the df filtered for J gene splices, changeo_H is
    the changeo_db from bracer summarise, and threshold is defines the ratio of
    max counts supporting the major J to C splice junction divided by the sum of all J to C splice junctions
    returns an SJout like df with only cells putatively in the act"""
    #get rid of cells expressing IGHD(the dual expression of M and D confounds my filters)
    _df = ab_sjoutdf[~ab_sjoutdf.exon_start.str.contains("IGHD")].copy()
    _df = _df[_df.exon_start.str.contains("exon1")]
    #divide the max unique mapping J to C splicing junction by the sum of all J to C splicing junctions
    _x = _df.groupby('cell').max()/_df.groupby('cell').sum()
    # cells in the act will have a number less than 1 for this value
    act_cells = _x[_x.unique_mapping < threshold].index
    #cells possibly in the act because they have multiple constant regions to J gene splicing events
    _act_cells_df = _df[_df.cell.isin(act_cells)]
    #cross reference to the assembly in order to use only the productive H chain transcript
    changeo_merge = changeodb_H[['CELL', 'J_CALL', 'IN_FRAME', 'FUNCTIONAL']]
    changeo_merge.columns = ['cell', 'J_CALL', 'IN_FRAME', 'FUNCTIONAL']
    
    merged_df = pd.merge(changeo_merge, _act_cells_df, on='cell', how='inner')
    productive_J_tx = merged_df[merged_df['J_CALL'].str.split('*', expand = True)[0] == merged_df['J_exon']]
    productive_J_tx = productive_J_tx[productive_J_tx.IN_FRAME == True]
    productive_J_tx = productive_J_tx[productive_J_tx.FUNCTIONAL == True]
    productive_J_tx = productive_J_tx[productive_J_tx.duplicated(subset=['cell', 'J_CALL'], keep=False)]
    #log transform counts
    productive_J_tx['log2_unique_mapping'] = np.log2(productive_J_tx['unique_mapping'])
    cells_in_act = productive_J_tx.drop_duplicates()
    #Plot
    if plot == True:
        f, ax = plt.subplots(figsize=(7, 20))

        ax = sns.barplot(data = cells_in_act, y = 'cell', x = 'log2_unique_mapping', hue = 'exon_start')
        plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad = 0)
    print(str(cells_in_act.cell.unique().shape[0]), 'cells in act of switching')
    return cells_in_act


def plotSwitchAndAbtx(ab_tx, switch_tx, cell_list):
    """makes a point plot of the IGH locus where each observation is a cell and an observation consists
    of counts for each of the genes at the IgH locus"""
    #Make ab_tx dataframe long form (cell, exon_start, unique_mapping)
    _dfab = ab_tx[ab_tx.exon_start.str.contains('exon1')].copy()
    _dfab_sum = _dfab.groupby(['cell', 'exon_start']).sum()
    _dfab_sum['unique_mapping_log2'] = np.log2(_dfab_sum['unique_mapping'])
    _dfab_long = _dfab_sum.unique_mapping_log2.unstack().fillna(0)


    new_columns = []
    for i in _dfab_long.columns:
        gene_name = i.split('_')[0]
        new_columns.append(gene_name)

    _dfab_long.columns = new_columns

    point_plot_df = _dfab_long

    point_plot_df = point_plot_df.reset_index()
    point_plot_df = pd.melt(point_plot_df, value_vars = point_plot_df.columns[1:], id_vars = point_plot_df.columns[0])
    point_plot_df.columns = ['cell', 'exon', 'log2CPM']
    point_plot_df['switch'] = 'Antibody'
    
    _dfab = point_plot_df
    
    
    ###
    #Switch df
    ####
    
    _dfswitch = switch_tx[switch_tx.exon_start.str.contains('exon1')].copy()
    _dfswitch_sum = _dfswitch.groupby(['cell', 'exon_start']).sum()
    _dfswitch_sum['unique_mapping_log2'] = np.log2(_dfswitch_sum['unique_mapping'])
    _dfswitch_long = _dfswitch_sum.unique_mapping_log2.unstack().fillna(0)
    new_columns = []
    for i in _dfswitch_long.columns:
        gene_name = i.split('_')[0]
        new_columns.append(gene_name)

    _dfswitch_long.columns = new_columns
    
    
    point_plot_df = _dfswitch_long

    point_plot_df = point_plot_df.reset_index()
    
    point_plot_df = pd.melt(point_plot_df, value_vars = point_plot_df.columns[1:], id_vars = point_plot_df.columns[0])
    point_plot_df.columns = ['cell', 'exon', 'log2CPM']
    point_plot_df['switch'] = 'Switch'
        
    _dfswitch = point_plot_df
    
    
    df = pd.concat([_dfswitch, _dfab])
    
    point_plot_df = df[df.cell.isin(cell_list)]

    ## Point Plot individual clones 
    #Data munging 
    
    #plotting 
    fig, ax = plt.subplots(1, 1, figsize=(12,4))

    sns.pointplot(data = point_plot_df,
                  x = 'exon', y = 'log2CPM', hue='switch',  dodge =.2, 
                  join = True, 
                  units = 'cell',
                  estimator = None,
                  #order by actual locus order
                  order = ['IGHM', 'IGHD', 'IGHG3','IGHG1', 'IGHA1', 
                           'IGHG2', 'IGHG4', 'IGHE', 'IGHA2']) 
    sns.despine()
    ax.set_ylabel("Log2_Unique_Junction_Reads")
    ax.set_xlabel("Genes")
    sns.set_context("talk")
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    sns.set_style("whitegrid", {'axes.grid' : False})
    
    fig, ax = plt.subplots(1, 1, figsize=(12,4))

    #sns.pointplot(data = IGH_locus_df,
                  #x = 'exon', y = 'log2CPM', hue = 'cell', dodge =.2, 
                  #join = True,
                  #order =    
    sinaplot(data= point_plot_df, x='exon', y='log2CPM',                   #order by actual locus order
                  order = ['IGHM', 'IGHD', 'IGHG3','IGHG1', 'IGHA1', 
                           'IGHG2', 'IGHG4', 'IGHE', 'IGHA2'])
    sns.despine()
    ax.set_ylabel("Log2_Unique_Junction_Reads")
    ax.set_xlabel("Genes")
    sns.set_context("talk")
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    sns.set_style("whitegrid", {'axes.grid' : False})

    fig, ax = plt.subplots(1, 1, figsize=(12,4))

    sns.pointplot(data = point_plot_df,
                  x = 'exon', y = 'log2CPM', dodge =.2, 
                  join = True,                   #order by actual locus order
                  order = ['IGHM', 'IGHD', 'IGHG3','IGHG1', 'IGHA1', 
                           'IGHG2', 'IGHG4', 'IGHE', 'IGHA2']) 
                      
    sns.despine()
    ax.set_ylabel("Log2_Unique_Junction_Reads")
    ax.set_xlabel("Genes")
    sns.set_context("talk")
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    sns.set_style("whitegrid", {'axes.grid' : False})
