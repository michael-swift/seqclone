import pandas as pd
import numpy as np
import scanpy as sc

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
