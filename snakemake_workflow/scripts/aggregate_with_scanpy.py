import scanpy as sc
import anndata as ad
# get files from snake make
files = snakemake.input
# params from snakemake
min_genes = int(snakemake.params.min_genes)
min_counts = int(snakemake.params.min_counts)

def load_and_filter(filename, min_genes=min_genes, min_counts=min_counts):
    adata = sc.read_h5ad(filename)
    # parse file name for sample_uid
    sample_uid=filename.split("/")[-1]
    sample_uid=sample_uid.split(".")[0]
    print(sample_uid)
    adata.obs["parsed_sample_uid"] = sample_uid
    adata.obs_names_make_unique()
    adata.var_names_make_unique()
    # filter very leniently
    sc.pp.filter_cells(adata, min_genes = min_genes)
    sc.pp.filter_cells(adata, min_counts = min_counts)
    return adata

adata_list = []

# construct aggregated object
for i, file_name in enumerate(files):
    adata_list.append(load_and_filter(file_name))
    print(i, " anndatas concatenated")

adata = ad.concat(adata_list)

print(str(snakemake.output))
adata.write_h5ad(str(snakemake.output))
print("Done!!")
