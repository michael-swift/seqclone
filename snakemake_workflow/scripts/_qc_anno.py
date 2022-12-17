import pandas as pd
import numpy as np
import scanpy as sc
get_ipython().run_line_magic('matplotlib', 'inline')
from matplotlib import pyplot as plt
import matplotlib as mpl
import seaborn as sns
import pathlib
import celltypist
from celltypist import models
params = {
    'font.size': 12,
    'axes.titlesize': 12,
    'axes.labelsize': 12,
    'legend.fontsize': 12,
    'xtick.labelsize': 8,
    'ytick.labelsize': 10,
    'font.family': "Arial",
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
    'figure.dpi': 100
   }

mpl.rcParams.update(params)

sns.set_style("ticks")
sns.set_context(context='paper')
savefig_args = {"dpi": 300, "bbox_inches": "tight", "pad_inches": 0, "transparent": True}
mpl.rc('savefig', dpi=300)
output_dir='figures/QCandAnnotation'
pathlib.Path(output_dir).mkdir(parents=True, exist_ok=True)
output_suffix = ""
output_formats = [".png", ".svg"]
sc.settings.figdir = output_dir
sc.set_figure_params(format='pdf', transparent=True,)
def save_figure(fig, name, output_dir=output_dir, output_suffix=output_suffix, output_formats=output_formats, savefig_args=savefig_args):
    for output_format in output_formats:
        fig.savefig(output_dir + "/" + name + output_suffix + output_format, **savefig_args)
    return None

pd.set_option('display.max_rows', 50)
pd.set_option('display.max_columns', 20) 
pd.set_option('display.width', 100)
get_ipython().run_line_magic('autoreload', '2')


def recluster(adata, batch_correct, sample_uid = 'sample_uid'):
    sc.pp.pca(adata)
    if batch_correct == True:
        sc.external.pp.bbknn(adata, batch_key=sample_uid)
    else:
        sc.pp.neighbors(adata, n_neighbors=10)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=0.2)
    sc.pl.umap(adata, color = sample_uid)
    return adata

import '../../scripts/plotting_helper.py')

adata = sc.read_h5ad('/home/michaelswift/repos/shared_data/BCD2/scirpy_processed_cellbender.h5ad')

adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=True, inplace=True)
sc.pl.violin(adata, ['n_genes_by_counts', 'log1p_total_counts', 'pct_counts_mt'], groupby='sample_uid', stripplot=False,
              multi_panel=True, rotation = 90, save = "_qc_pre_filter")

# filter cell
sc.pp.filter_cells(adata, min_genes=500)
sc.pp.filter_cells(adata, min_counts=2000)
sc.pp.filter_cells(adata, max_counts=100000)
# filter mt 
adata = adata[adata.obs.pct_counts_mt < 10]
sc.pl.violin(adata, ['n_genes_by_counts', 'log1p_total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True)
sc.pl.violin(adata, ['n_genes_by_counts', 'log1p_total_counts', 'pct_counts_mt', "n_genes"], groupby='sample_uid', stripplot=False,
              multi_panel=True, rotation = 90, save = "_quality_metrics_post_filter")

adata.layers['counts'] = adata.X

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata, base = 2, chunk_size=10000 )
sc.pp.highly_variable_genes(adata, n_top_genes=4000)
sc.pl.highly_variable_genes(adata)
adata.raw = adata
sc.pp.scale(adata, max_value = 10)
# remove IGH and IGL variable genes from highly variable genes for clustering analysis 
adata.var.loc[adata.var.index.str.contains("IGHV|IGLV|IGKV"), 'highly_variable'] = False

adata = recluster(adata, batch_correct=True)
adata.obs_names_make_unique()
# rename lane for interpretability:
rename_dict = {'10K_BMMNC' : 'BMMNC', '10K_PBMC': 'PBMC',
               '20K_PBMC': 'PBMC' , "BMPC_L1" : "BM CD138+", "BMPC_L2": "BM CD138+",
               'L01_Input': 'Day 0', 'L04_T4': "Day 4", "L06_T8":"Day 8", "L08_T12" : "Day 12"}   

adata.obs['sample_id'] = adata.obs.sample_uid.map(rename_dict)
anno_1 = "celltypist_all"
#Download all the available models.
models.download_models()
#Provide the input as an `AnnData`.
predictions = celltypist.annotate(adata, model = 'Immune_All_Low.pkl', majority_voting=True)
# Celltypist Annotations to bcells object
adata.obs[anno_1] = predictions.predicted_labels.majority_voting

sc.pl.umap(adata, color=anno_1, frameon=True, 
          title="Celltypist Annotation", add_outline=True, legend_fontoutline=2, legend_fontsize='xx-small', save="_CellType_Annotation")

sc.pl.umap(adata, color="sample_uid", frameon=True,
          title="Sample annotation", add_outline=True, legend_fontoutline=4, legend_fontsize='x-small', save="_sample_uid_Annotation")

sc.pl.umap(adata, color="CD1C", frameon=False,
          title="Sample Annotation", add_outline=True, legend_fontoutline=4, legend_fontsize='x-small', save="CD1C")

# remove NK Cells which become contaminants downstream and appear to confuse Celltypist
sc.tl.score_genes(adata, gene_list=['NCAM1',
 'GZMB',
 'VIM',
 'NKG7',
 'PRF1',
 'RNF213',
 'GNLY',
 'CTSW',
 'CST7',
 'GZMA'])

sns.displot(adata.obs.score)
sc.pl.umap(adata, color = 'score')

adata = adata[adata.obs.score < 0.1]


# In[21]:


sc.pl.umap(adata, color="sample_uid", frameon=False,
          title="Sample annotation", add_outline=True, legend_fontoutline=4, legend_fontsize='x-small', save="adata_sample_uid_Annotation")


# # Annotate The B Cells with Subtypes

# In[22]:


bcells = adata[adata.obs[anno_1].str.contains('B cells|Plasma')]
bcells.obs.sample_uid.value_counts()


# In[23]:


# Make all IG not a highly variable gene for clustering purposes
bcells.var.loc[bcells.var.index.str.contains("IGH|IGL|IGK"), 'highly_variable'] = False


# In[24]:


bcells = recluster(bcells, batch_correct=True)


# In[25]:


# Celltypist Annotations to bcells object
anno_2 = 'celltypist'
predictions = celltypist.annotate(bcells, model = 'Immune_All_Low.pkl', majority_voting=True)
bcells.obs[anno_2] = predictions.predicted_labels.majority_voting


# In[26]:


sc.tl.leiden(bcells, resolution=1)
sc.pl.umap(bcells, color = 'leiden')



# In[27]:


sc.pl.umap(bcells, size= 10, color = anno_1, edges=True, edges_width=0.01)


# ## Add columns to object for use in creating figures

# In[28]:


bcells = bcells[bcells.obs[anno_2].str.contains('B cells|Plasma')]


# In[29]:


bcells.obs['isotype_simple'] = bcells.obs.IR_VDJ_1_c_call.map(IGH_simplify())


# In[30]:


bcells.obs['fraction_mutated_bases'] = 1-bcells.obs['IR_VDJ_1_v_identity']


# In[ ]:


# filtering out low quality B cells
# Filters out some probably artefactual cells
bcells = bcells[bcells.obs.has_ir == 'True']
# Filters out putatively dead / apoptotic B cells
# filter B cells without a Constant Region Call
bcells = bcells[bcells.obs.dropna(subset='IR_VDJ_1_c_call').index, :]

bcells = bcells[bcells.obs.dropna(subset='isotype_simple').index, :]

# filter out weird cells where they have a light chain called as the heavy chain
bcells = bcells[~bcells.obs.IR_VDJ_1_c_call.str.contains("IGK|IGL"),:]
# filter out B cells without mutational information
bcells = bcells[bcells.obs.simple_mutation_status.dropna().index]
#bcells.obs['fraction_mutated_bases'] = 1 - bcells.obs['IR_VDJ_1_v_identity']

bcells.obs['bcelltype'] = bcells.obs['celltypist']
cmap = plt.get_cmap('rainbow')
colors = cmap(np.linspace(0, 1, len(bcells.obs.bcelltype.cat.categories)))
bcells.uns["bcelltype_colors"] = list(map(mpl.colors.to_hex, colors))
bcells.uns["bcelltype_colors_dict"] = dict(zip(bcells.obs.bcelltype.unique().sort_values(), bcells.uns['bcelltype_colors']))
color = 'bcelltype'
palette = bcells.uns[color+'_colors']

bcells.uns['celltype_colors'] = None

bcells.obs['bcelltype'] = bcells.obs.celltypist.str.replace("Age-associated ", "").values

bcells.obs.bcelltype = bcells.obs['bcelltype'].str.replace("germinal center", "GC")
bcells.obs.bcelltype = bcells.obs['bcelltype'].str.replace("Proliferative", "Prolif.")

# combine mutation status and isotype call
bcells.obs['combined_ab'] = bcells.obs['mutation_status'].astype(str) + "_" + bcells.obs['isotype_simple'].astype(str)
bcells.obs['bcelltype_multi'] = bcells.obs['bcelltype'].astype(str) + "_" + bcells.obs['simple_mutation_status'].astype(str) + "_" + bcells.obs['isotype_simple'].astype(str)

bcells.obs['switched'] = bcells.obs.IR_VDJ_1_c_call.map(IGH_switched())


# # what's with all the ABCs?

# In[ ]:


models.models_description()

model = models.Model.load(model = '/home/michaelswift/.celltypist/data/models/Immune_All_Low.pkl')
top_genes = model.extract_top_markers("Age-associated B cells", 15)
sc.pl.umap(bcells, color = list(top_genes) + ['bcelltype'])
# these genes doen't appear specific, a B cell


# In[ ]:


sc.pl.umap(bcells, color = ['ITGAX', 'CD1C', 'CD200', "CXCR5", 'BCL6', 'MALAT1', "AICDA", "MYC", "EBI3", 'sample_id', 'mutation_status'], use_raw = False, size = 12)


# In[ ]:


sc.pl.umap(bcells, color = ['ITGAX', 'CD1C', 'CD200', "CXCR5", 'BCL6', 'MALAT1', "CD38", "CCR6",'sample_id', 'mutation_status'], use_raw = True, size = 10)


# In[ ]:


bcells.write_h5ad('../../processed_data/h5ad_objects/bcells.12.1.2022.h5ad')

