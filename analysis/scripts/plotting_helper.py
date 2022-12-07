### Dictionaries for color definition etc.
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests
from matplotlib import pyplot as plt
import seaborn as sns
import scanpy as sc
from matplotlib.lines import Line2D
###########################
## Plotting Dictionaries ##
###########################
def celltype_colors():
    color_dict = {'Plasma B cell': "#29975D" , 'Myeloid cell': "#972963",
                    'Resting B cell': "#6BB9CB",
       'Activated B cell': "#000000", 'T cell': "#EF9A9A", 'NK / Dead Cell':"#ABABAB" , 'NK cell': "#8A679A"}
    return color_dict

def sample_id_colors():
    
    sample_id_color_dict = {'BM CD138+': '#1f09ff',
             'BMMNC': '#0977ff',
             'Day 0': '#14e0be',
                            'Day 12':"#ff3939",
                            'Day 8':"#da5435",
                            'Day 4':"#c0cd30",
                            'PBMC':"#09a3ff"                         }
    return sample_id_color_dict

def switched_colors():
    color_dict = {"IGHM|IGHD": "#14e0be",
                  "switched":"#000000"}
    
    return color_dict

def mutation_colors():
    mutation_colors = {'mutated': '#FF5733', 'germline' : '#004D40', 'heavily mutated':'#581845'}
    return mutation_colors

def IGH_colors():
    
    IGH_color_dict = {'IGHA1': '#FFAB91',
             'IGHA2': '#FFD180',
             'IGHD': '#00E676',
             'IGHM': '#69F0AE',#"#00ff7f",
             'IGHG2': '#B388FF',
             'IGHG4': '#008784',
             'IGHG1': '#0B3954',
             'IGHG4': '#6A1B9A',
             'IGHG3': '#536DFE',
             'IGHE':'#FF5A5F'}
    return IGH_color_dict

def IGH_simplify():
    
    IGH_simple = {'IGHA1': 'IGHA',
             'IGHA2': 'IGHA',
             'IGHD': 'IGHD',
             'IGHM': 'IGHM',#"#00ff7f",
             'IGHG2': 'IGHG',
             'IGHG4': 'IGHG',
             'IGHG1': 'IGHG',
             'IGHG4': 'IGHG',
             'IGHE':'IGHE'}
    return IGH_simple

def IGH_colors_simple():
    
    IGH_color_dict = {'IGHA': '#e0b96a',
             'IGHD': '#798658',
             'IGHM': '#17a589',#"#00ff7f",
             'IGHG': '#0D47A1',
             'IGHE':'#8b4cf4'}
    return IGH_color_dict

def plot_order():
    order = ['BMMNC',
 'PBMC',
 'BM CD138+',
 'Day 0',
 'Day 4',
 'Day 8',
 'Day 12']
    return order

def timecourse():
    order = [
 'Day 0',
 'Day 4',
 'Day 8',
 'Day 12']
    return order

def pseudo_timecourse():
    order = [
 'Day 0',
 'Day 4',
 'Day 8',
 'Day 12', 'BM CD138+']
    return order

mutation_cutoffs = [0.02, 0.10]

def IGH_switched():
    
    IGH_switched = {'IGHA1': 'switched',
             'IGHA2': 'switched',
             'IGHD': 'IGHM|D',
             'IGHM': 'IGHM|D',#"#00ff7f",
             'IGHG2': 'switched',
             'IGHG4': 'switched',
             'IGHG1': 'switched',
             'IGHG4': 'switched',
             'IGHE':'switched'}
    return IGH_switched

def timepoint_shapes():
    timepoint_shapes_dict = {"Day 0": "triangle", "Day 4": "square", "Day 8": "pentagon", "Day 12": "octagon"}
    return timepoint_shapes_dict

def celltype_colors():
    color_dict = {'Plasma B cell': "#29975D" , 'Myeloid cell': "#972963",
                    'Resting B cell': "#6BB9CB",
       'Activated B cell': "#000000", 'T cell': "#EF9A9A", 'NK / Dead cell':"#ABABAB" , 'NK cell': "#8A679A"}
    return color_dict


def bcelltype_colors_dict():
    colors = {'B cells' : '#8000ff', 'Memory B cells' : '#1996f3', 'Naive B cells':'#4df3ce', 'Plasma cells':'#b2f396', 'Plasmablasts': '#ff964f', 'Prolif. GC B cells': '#ff0000'}
    return colors


IGH_order = {'IGHM':0, 'IGHD':1, 'IGHG3':2, 'IGHG1':3, 'IGHA1':4, "IGHG2":5, 'IGHG4':6, "IGHE":7, "IGHA2":8}

IGH_simple_markers = {'IGHM':'^', 'IGHD':"v", 'IGHG':"o", 'IGHA':"s","IGHE":"*"}
UMAP_variables_of_interest = ['sample_id', 'isotype_simple', 'switched', 'simple_mutation_status']

###########################
##        Execute        ##
###########################

timepoint_shapes = timepoint_shapes()
celltype_colors = celltype_colors()
bcelltype_colors = bcelltype_colors_dict()
sample_id_colors = sample_id_colors()
switched_colors = switched_colors()
mutation_colors = mutation_colors()
plot_order = plot_order()
timecourse = timecourse()
pseudo_timecourse = pseudo_timecourse()
igh_colors = IGH_colors()
igh_colors_simple = IGH_colors_simple()
igh_genes = list(IGH_order.keys())
# other settings
savefig_args = {"dpi": 300, "bbox_inches": "tight", "pad_inches": 0, "transparent": True}
output_suffix = ""
output_formats = [".png", ".svg"]

###########################
##     Helper Fxns       ##
###########################
def selection_helper(df, param, value):
    selector = df[param].value_counts() > value
    idxs = selector[selector == True].index
    return(df[df[param].isin(idxs)])

def shuffleWithinGroup(data, label, group_name):
    """ returns a shuffled series of the label (column name), shuffled within the groups supplied """
    list_of_dfs = []
    for group, frame in data.groupby(group_name):
        frame.loc[:,label] = np.random.permutation(frame.loc[:,label])
        list_of_dfs.append(frame)
    df = pd.concat(list_of_dfs)
    return df.loc[:,label]

def correct_p(df, column, method):
    result = multipletests(df[column], method = method)
    df['significant'] = result[0]
    df[column + 'corrected_pvalue'] = result[1]
    return df

def bootstrap_routine(df, group, label, boots, samplesize):
    original = df[[group, label]]
    print(group)
    observed = pd.DataFrame(original.groupby(group)[label].value_counts(normalize = True, dropna = False))
    observed.reset_index(inplace = True)
    observed.columns = [group, label, 'proportion']
    # take all samples -> then split into distinct bootstraps
    print(boots, samplesize)
    bootstrapped_samples = original.groupby(group).sample(n = boots * samplesize, replace = True)
    # reshape
    data = bootstrapped_samples.values.reshape(-1, boots)
    
    df = pd.DataFrame(data)
    bootstrapped_list_of_stat_df = []
    # index bootstrapped samples to create single bootstrapped dataframe
    for i in np.arange(df.shape[1], step = 2):
        bootstrap = pd.DataFrame([df.iloc[:,i], df.iloc[:,i+1]]).T
        bootstrap.columns = [group, label]
        # perform statistical calculation on bootstrapped dataframe
        dfs = []
        for index, grp in bootstrap.groupby(group):
            _df = pd.DataFrame((bootstrap.groupby(group)[label].value_counts(normalize = True).xs('Day 0') / bootstrap.groupby(group)[label].value_counts(normalize = True).xs(index)) ** -1 )
            _df['sample_id'] = index
            dfs.append(_df)

        data = pd.concat(dfs)
        data.dropna(subset = ['mutation_status'], inplace = True)
        #data.reset_index(inplace = True)
        bootstrapped_list_of_stat_df.append(data)
    df = pd.concat(bootstrapped_list_of_stat_df)
    df.columns = ['ratio', 'sample_id']
    df.reset_index(inplace = True)
    return df, observed

# visualize the color dictionaries
def plot_colortable(colors, title, sort_colors=True, emptycols=0):

    cell_width = 230
    cell_height = 20
    swatch_width = 48
    margin = 20
    topmargin = 40

    # Sort colors by hue, saturation, value and name.
    by_hsv = [(v, k) for k, v in colors.items()]
    
    if sort_colors is True:
        by_hsv = sorted(by_hsv)
    names = [name for hsv, name in by_hsv]

    n = len(names)
    ncols = 2 - emptycols
    nrows = n // ncols + int(n % ncols > 0)

    width = cell_width * 4 + 2 * margin
    height = cell_height * nrows + margin + topmargin
    dpi = 72

    fig, ax = plt.subplots(figsize=(width / dpi, height / dpi), dpi=dpi)
    fig.subplots_adjust(margin/width, margin/height,
                        (width-margin)/width, (height-topmargin)/height)
    ax.set_xlim(0, cell_width * 4)
    ax.set_ylim(cell_height * (nrows-0.5), -cell_height/2.)
    ax.yaxis.set_visible(False)
    ax.xaxis.set_visible(False)
    ax.set_axis_off()
    ax.set_title(title, fontsize=24, loc="left", pad=10)

    for i, name in enumerate(names):
        row = i % nrows
        col = i // nrows
        y = row * cell_height

        swatch_start_x = cell_width * col
        swatch_end_x = cell_width * col + swatch_width
        text_pos_x = cell_width * col + swatch_width + 7

        ax.text(text_pos_x, y, name, fontsize=12,
                horizontalalignment='left',
                verticalalignment='center')

        ax.hlines(y, swatch_start_x, swatch_end_x,
                  color=colors[name], linewidth=18)

    return fig
def plot_cascade(adata, gene, label, hue, ax, shuffle=False, clone_size=4, palette=sample_id_colors, size = 3):
    df = sc.get.obs_df(adata, keys=[gene], use_raw=True)
    data = adata
    # Only Clones
    if label == 'clone_id':
        _selector = data.obs[label].value_counts() > clone_size
        _selector = _selector.index[_selector == True]
        data = data[data.obs[label].isin(_selector)]

    data.obs.loc[:, gene] = df[gene]
    if shuffle == True:
        data.obs[label] = shuffleWithinGroup(data.obs, label, hue)
    data = pd.DataFrame([data.obs[label].to_list(), data.obs[gene].to_list(), data.obs[hue].to_list()])
    data = data.T
    data.columns = [label, gene, hue]
    # calculate average CV:
    avg_CV = np.round((data.groupby(label)[gene].std() / data.groupby(label)[gene].mean()).mean(), decimals=2)
    print("{} \n average CV: {}".format("Shuffled: " + str(shuffle),  avg_CV))
    ax = sns.stripplot(data = data, x = label, y = gene, size = size, alpha = 1, linewidth=0.5, 
                  order = data.groupby(label)[gene].mean().sort_values().index, hue = hue, palette = palette, ax = ax)
    #ax.text(1, (data[gene].max() - data[gene].min())/ 1.25, 'Avg CV: {}'.format(avg_CV), style='italic', bbox={
    #    'facecolor': 'white', 'alpha': 0.8, 'pad': 1})
    
    return ax

def plot_dual_cascade(adata, label, hue, palette, clone_size, gene, width, height):
    shuffled = [False, True]
    fig = plt.figure(figsize = (width*2, height/2))
    gs = fig.add_gridspec(1,2, hspace=0)
    axs = gs.subplots(sharex=True, sharey=True)
    for d, ax in enumerate(axs):
        plot_cascade(adata, gene, label, hue, ax=ax, shuffle=shuffled[d], clone_size=clone_size, palette=palette)
        ax.label_outer()
        ax.set_xlabel('clones')
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0)
        # for major ticks
        ax.set_xticks([])
        # for minor ticks
        #ax.set_xticks([], minor=True)
        if d == 0:
            ax.get_legend().remove()
            ax.set_ylabel("{}\n(log umi\n per 10K)".format(gene))
        if shuffled[d] == True:
            ax.set_title('Permuted Labels')
        if shuffled[d] == False:
            ax.set_title('True Labels')
    #save_figure(fig, "{}_cascade".format(gene))
    return fig

def _prepare_seqclone_df(df):
    # filter out genes which seem unreasonable to analyze:
    df = df[df.n_cells_by_counts > 10]
    df = df.dropna(subset = 'pvalue')
    df = df[~df.index.str.contains('AC1|AP0|AC2|AC0')]
    df = df[df.log1p_mean_counts > 0.005]
    df['normalized_effect_size'] = df['effect_size'] / df['mean_counts']
    df.sort_values(by='normalized_effect_size', inplace= True)
    df.rename(columns={'freq_pvalue':'nulls'}, inplace = True)
    # p value plus a pseudocount
    df['frequentist_p_value'] = (df.nulls + 1) / 10000
    df = _correct_p(df, 'frequentist_p_value', method='fdr_bh')
    #df['significant'] = df.frequentist_p_value < 0.01
    df['absolute_effect_size'] = df.effect_size.abs()
    df['normalized_absolute_effect_size'] = df.absolute_effect_size / df.log1p_mean_counts
    df['known clonal gene'] = df.index.str.contains('IGKV|IGHV|IGLV')
    df['significant-IgVariable'] = df.significant.astype(str) + '-' + df.index.str.contains('IGKV|IGHV|IGLV').astype(str)
    df['legend'] = df['significant-IgVariable'].map({'True-True':"Known Clonal & p < 0.01", "True-False": "p < 0.01", "False-False": "p > 0.01"})
    return df

def _correct_p(df, column, method):
    result = multipletests(df[column], method = method, alpha = 0.01)
    df['significant'] = result[0]
    df[column + '_corrected'] = result[1]
    return df

def plot_volcano(df, pvalue_threshold, effect_cutoff, genes_to_label, x, y, xplacement, yplacement):
    g = sns.jointplot(data = df, x = x , y = y, hue='legend', 
                      palette='Accent', joint_kws={"s":12, "edgecolor":"none", 'alpha':0.5}, marginal_kws={'log_scale':True, "cut":0, }, marginal_ticks=True, rasterized = True, ratio = 3)
    g.set_axis_labels(xlabel = 'q-value', ylabel = 'Clonal Index')
    g.fig.axes[0].hlines(xmax=1.1, xmin=0, y = effect_cutoff, color = 'k', linestyle = '-')
    g.fig.axes[0].vlines(ymax=1, ymin=0, x = pvalue_threshold, color = 'k', linestyle = '-')
    # make equal bins
    genes = df.loc[genes_to_label,:]
    for i, gene in enumerate(genes.index):
        x_co = genes.loc[gene,x]
        y_co = genes.loc[gene,y]
        g.fig.axes[0].annotate(gene, (x_co, y_co), xytext=(xplacement[i] , yplacement[i]),
        bbox=dict(boxstyle="round", alpha=0.1), 
        arrowprops = dict(arrowstyle="simple", mutation_scale=0.1))
    return g


def getEmbeddingCoordinates(adata, embedding):
    embedding_key = 'X_' + embedding
    _df_embedding = sc.get.obs_df(adata, obsm_keys = [(embedding_key, 0), (embedding_key, 1)])
    return _df_embedding

def plot_ClonalUMAP(adata, embedding, clone_label, size_cutoff, point_size, size_multiple, custom_legend):
    pre_dict = adata.obs.clone_id.value_counts().reset_index().iloc[:,:1]

    mapper = pre_dict.to_dict()['index']

    inv_map = {v: k for k, v in mapper.items()}

    adata.obs['CloneLabel'] = adata.obs[clone_label].map(inv_map)

    selector = adata.obs.CloneLabel.value_counts()
    selector = selector[selector > 1].index
    adata.obs['IN_CLONE'] = adata.obs.CloneLabel.isin(selector)

    selector = adata.obs.CloneLabel.value_counts()

    selector = selector[selector > size_cutoff].index
    adata.obs['IN_BIG_CLONE'] = adata.obs.CloneLabel.isin(selector) 

    df = getEmbeddingCoordinates(adata, 'umap')

    df['CloneLabel'] = adata.obs.CloneLabel
    df['IN_CLONE'] = adata.obs.IN_CLONE
    df['IN_BIG_CLONE'] = adata.obs.IN_BIG_CLONE
    
    sns.set_palette('Accent')
    clone_color = 'lightcoral'
    non_clone_color = 'k'
    fig, ax = plt.subplots(1,1, figsize = (2.5,2.5)) 
    s = point_size
    # Clones
    clones = df[df.IN_CLONE == True]
    # Big Clones
    BigClones = df[df.IN_BIG_CLONE == True]
    # Hacky way to reduce number of categories which seems to be inherited from the superset df
    BigClones = pd.DataFrame(BigClones.values)
    BigClones.columns = clones.columns
    BigClones['CloneLabel'] = BigClones.CloneLabel.astype('category')
    BigClones['Clone Label'] = 'Clone ' + BigClones.CloneLabel.astype(str)
    # Non Clones
    non_clones = df[df.IN_CLONE == False]
    x = non_clones.iloc[:,0]
    y = non_clones.iloc[:,1]
    ax = sns.scatterplot(x = x , y = y, color = 'k', alpha=0.6, s=s, rasterized = True)
    # Clones
    x = clones.iloc[:,0]
    y = clones.iloc[:,1]
    ax = sns.scatterplot(x = x , y = y, color = 'lightcoral', edgecolor = 'k', linewidths = 0.5, alpha=0.01, s=s*0.8, rasterized = True)
    # Big Clones
    x = BigClones.iloc[:,0]
    y = BigClones.iloc[:,1]
    ax = sns.scatterplot(data = BigClones, x = x, y = y, alpha = 1, hue='Clone Label', linewidths = 0.5, edgecolor="k", s = s*size_multiple, palette='Dark2', rasterized = True)
    if custom_legend:
        legend_elements = [Line2D([0], [0], marker='o', color=clone_color, label='Clonal Family Member\n Detected',
                          markerfacecolor=clone_color, linewidth = 0, markersize=4,), 
                   Line2D([0], [0], marker='o', linewidth = 0, markersize=4, color=non_clone_color,label='Clonal Family Member\n not Detected',
                          markerfacecolor=non_clone_color)]
        ax.legend(handles=legend_elements, bbox_to_anchor=(1.1, 1), loc='upper left', borderaxespad=0)
    else:
        plt.legend(loc = 5, bbox_to_anchor = (1.5, 1))
    # Remove UMAP units labels
    ax.axes.xaxis.set_ticklabels([])
    ax.axes.yaxis.set_ticklabels([])

    plt.xlabel('UMAP 1')
    plt.ylabel('UMAP 2')
    return fig