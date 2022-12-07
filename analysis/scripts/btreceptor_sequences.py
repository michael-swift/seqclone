from __future__ import division
import pandas as pd
import numpy as np
import Levenshtein
from scipy.spatial.distance import squareform
from scipy.sparse.csgraph import connected_components
from itertools import combinations

def df_pw_edit(frame):
    """ Returns array of pairwise edit distances in square form """

    ed = np.zeros(int((frame.shape[0]/2)*(frame.shape[0]-1)), dtype='float')
    for c, (x, y) in enumerate(combinations(frame.cdr3aa.values, 2)):
        ed[c] = Levenshtein.distance(x, y) / np.max([len(x), len(y)])
    sq = squareform(ed)

    return sq


def df_lineages_from_subset(frame, similarity_cutoff):
    """ Returns an array of lineage membership based on a CDR3 cutoff """

    edit_sq = df_pw_edit(frame)

    n_groups, labels = connected_components(edit_sq <= round(1 - similarity_cutoff, 4))

    return labels


def df_add_lineages(dfin, similarity_cutoff):
    """ Returns input dataframe with additional lineage column
    Args:
        similarity_cutoff (float): e.g. 0.8 for 80% minimum cdr3aa similarity
    """

    dfin = dfin.copy()

    # unique index required for join
    if not dfin.index.is_unique:
        print("Input DataFrame index not unique, applying reset_index().")
        dfin.reset_index(drop=True, inplace=True)

    lincnt = 0
    lins = []

    for (v, j, _), sub in dfin.groupby(['v_call_no_allele',
                                       'j_call_no_allele',
                                       'cdr3aa_len']):
        if sub.shape[0] > 1:
            # CDR3 distance comparison
            sub_lineages = df_lineages_from_subset(sub, similarity_cutoff)
            lins += zip(sub.index, sub_lineages + lincnt)
            lincnt += np.unique(sub_lineages).shape[0]
        else:
            # single sequence belongs in its own lineage
            lins.append((sub.index.values[0], lincnt))
            lincnt += 1

    # adds a "lineage" column corresponding to the lineage number for that cell
    lins = pd.DataFrame(lins, columns=['index', 'lineage']).set_index('index')
    if 'lineage' in dfin.columns:
        dfin = dfin.drop('lineage', axis=1).join(lins)
    else:
        dfin = dfin.join(lins)

    return dfin

# Add lineages with in house pipeline to test
def adata_add_lineages_from_IR(adata, similarity_cutoff):
    
    df = adata.obs

    df['v_call_no_allele'] = df['IR_VDJ_1_v_call'].str.split("*", expand = True)[0]
    df['j_call_no_allele'] = df['IR_VDJ_1_j_call'].str.split("*", expand = True)[0]

    df['cdr3aa_len'] = df['IR_VDJ_1_junction_aa'].str.len()
    df['cdr3aa'] = df['IR_VDJ_1_junction_aa']

    df = df_add_lineages(df, similarity_cutoff=similarity_cutoff)

    adata.obs = df
    return adata

def df_shuffle_within_group(ird, label, group_name):
    """ returns a shuffled series of the label, shuffled within the groups supplied """
    list_of_dfs = []
    for group, frame in ird.groupby(group_name):
        frame.loc[:,label] = np.random.permutation(frame.loc[:,label])
        list_of_dfs.append(frame)
    df = pd.concat(list_of_dfs)
    return df.loc[:,label]

def selection_helper(df, param, value):
    selector = df[param].value_counts() > value
    idxs = selector[selector == True].index
    return(df[df[param].isin(idxs)])
