import scirpy as ir
import scanpy as sc

h5ad = str(snakemake.input)

# read adata
adata = sc.read_h5ad(h5ad)
# scirpy defining clonotypes
adata.obs_names_make_unique()
ir.tl.chain_qc(adata)
ir.pp.ir_dist(adata, n_jobs=20)
ir.tl.define_clonotypes(adata, receptor_arms="all", dual_ir="primary_only")

adata.write_h5ad(str(snakemake.output))
