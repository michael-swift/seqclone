### BtReceptor Edited for T cell data 

from __future__ import division
import pandas as pd
import numpy as np
import Levenshtein
from matplotlib import pyplot as plt
import matplotlib as mpl
from scipy.spatial.distance import squareform
from scipy.sparse.csgraph import connected_components
from itertools import combinations


from numpy.random import seed


def df_generate_node_dict(frame, singletons=False):
    """ Generates a node property dictionary for downstream graph-tool plotting
        Args:
            frame (pd.Dataframe) with the following columns:
                - clone_id (mandatory)
                - node_color (optional)
                - node_shape (optional)
                - node_size (optional)
                - node_stroke (optional)
            singletons (bool): whether to include cells that do not form
                multi-member clonal families
        Returns:
            dict: nested as follows: {node_id: {'property1': 'value1', ...}}
    """
    frame = frame.copy()

    # use some defaults if property columns are missing
    req_key_defaults = {'node_color': 'blue',
                        'node_shape': 'circle',
                        'node_size': 1,
                        'node_stroke': 0.5}
    for key, default in req_key_defaults.items():
        if key not in frame.columns:
            frame[key] = default  # add column to dataframe

    node_props = {}
    node_id = 0
    for lin, g in frame.groupby('clone_id'):
        if g.shape[0] > 1 or singletons:
            if g.shape[0] > 1:
                # add a germline to root all cells in a clone_id
                gid = node_id
                node_props[gid] = {'color': 'k',
                                   'ancestor': None,
                                   'size': 0.1,
                                   'shape': 'circle',
                                   'stroke': 0}
                node_id += 1
            else:
                # no germline root for a singleton
                gid = None

            # add cell(s) with properties
            for _, row in g.iterrows():
                node_props[node_id] = {'ancestor': gid,
                                       'color': row['node_color'],
                                       'shape': row['node_shape'],
                                       'size': row['node_size'],
                                       'stroke': row['node_stroke']
                                       }
                node_id += 1

    return node_props


def draw_gviz(node_dict, size_multiple=1, random_seed=42, **kwargs):
    """ Draw clonal network using graph-tool
    More information: graphtool edge / vertex parameters and examples:
        https://graph-tool.skewed.de/static/doc/draw.html#graph_tool.draw.graph_draw
        http://ryancompton.net/2014/10/05/graph-tools-visualization-is-pretty-good/
    Args:
        node_dict (dict): nested dictionary of node properties
            Generate this using df_generate_node_dict()
        size_multiple (int): scaling factor for node size (for convenience)
        **kwargs: keyword arguments passed to gt.graph-draw()
            e.g. output='file.pdf', layout='neato', output_size=(300,300)
    """
    import graph_tool.all as gt

    g = gt.Graph()
    vsizes = g.new_vertex_property("int")
    vcolors = g.new_vertex_property('string')
    vshapes = g.new_vertex_property('string')
    vpenwidth = g.new_vertex_property("float")  # stroke

    for node_id, node_props in node_dict.items():
        g.add_vertex()

        vshapes[g.vertex(node_id)] = node_props['shape']
        vcolors[g.vertex(node_id)] = node_props['color']
        vsizes[g.vertex(node_id)] = node_props['size']*size_multiple
        vpenwidth[g.vertex(node_id)] = node_props['stroke']

        # add edge to ancestor
        if node_props['ancestor'] is not None:
            g.add_edge(node_props['ancestor'], node_id)

    # seeds enable graph reproduction
    seed(random_seed)
    gt.seed_rng(random_seed)

    gt.graph_draw(g,
                  vertex_size=vsizes,
                  vertex_fill_color=vcolors,
                  vertex_shape=vshapes,
                  vertex_pen_width=vpenwidth,
                  vertex_color='k',  # stroke color
                  bg_color=[1, 1, 1, 1],  # white
                  edge_end_marker='none',
                  **kwargs)


def _detect_dataframe_cell_type(df):
    """ Detects whether cell type being analyzed is B or T cells
    """

    if df.v_call.str.contains('IG\wV').all():
        return 'B cell'
    elif df.v_call.str.contains('TR\wV').all():
        return 'T cell'
    else:
        raise ValueError('Mixed / unknown cell types found in dataframe: '
                         '{}'.format(df.v_call.str.slice(0, 2).unique()))


def rename_changeo_columns(df, revert=False):
    """ Simplifies dataframe column names.
    Args:
        df (pd.DataFrame)
        reverse (bool): whether to revert column names to original
    Returns:
        Dataframe
    """

    df = df.copy()

    rename_dict = {'sequence_id': 'seqid',
                   'cdr3_igblast_nt': 'cdr3nt',
                   'cdr3_igblast_aa': 'cdr3aa'}
    rename_dict.update({'{}_seq_length'.format(x):
                        '{}_len'.format(x) for x in ['v', 'd', 'j']})
    rename_dict.update({'{}_seq_start'.format(x):
                        '{}_start'.format(x) for x in ['v', 'd', 'j']})

    if not revert:
        # rename columns for convenience
        df.rename(columns=lambda x: x.lower(), inplace=True)
        df.rename(columns=rename_dict, inplace=True)
    else:
        # revert to original columns
        df.rename(columns={v: k for k, v in rename_dict}, inplace=True)
        df.rename(columns=lambda x: x.upper(), inplace=True)

    if 'cdr3aa' not in df.columns:
        print('\nWarning: change-o MakeDB.py was not run with "--cdr3".'
              ' Translating the amino acid CDR3 from IMGT nucleotide'
              ' sequence instead.')
        df['cdr3aa'] = df.cdr3_imgt.apply(
            lambda x: x if pd.isnull(x) else str(Seq(x).translate()))

    return df


def load_changeo_igblast_makedb(infile):
    """ Loads tab separated Change-O output into pandas dataframe. Adds
        additional convenience columns e.g. for clustering
    """

    ig = pd.read_csv(infile, sep='\t')

    ig = rename_changeo_columns(ig)

    cell_type = _detect_dataframe_cell_type(ig)
    if cell_type == 'B cell':
        ig['heavy'] = ig.v_call.str.startswith('IGH')
    elif cell_type == 'T cell':
        ig['heavy'] = (ig.v_call.str.startswith('TRB') |
                       ig.v_call.str.startswith('TRD'))

    # update / add columns
    ig['j_end'] = ig.j_start + ig.j_len
    ig['j_start_vdj'] = ig.j_start - ig.v_start
    ig['v_call_no_allele'] = ig.v_call.str.split('*', expand=True)[0]
    ig['j_call_no_allele'] = ig.j_call.str.split('*', expand=True)[0]
    ig['cdr3aa_len'] = ig.cdr3aa.str.len()

    # cast columns as bool if not already
    bool_cols = ['functional', 'indels', 'stop', 'in_frame']
    for col in bool_cols:
        if not pd.api.types.is_bool_dtype(ig[col]):
            ig[col] = ig[col] == 'T'

    return ig