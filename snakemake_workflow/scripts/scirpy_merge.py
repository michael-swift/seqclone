import pandas as pd
import scanpy as sc
import scirpy as ir

# I/O
gex_file = snakemake.input['h5']
airr_imm = snakemake.input['airr_imm']
json_TXG = snakemake.input['json_txg']
merged_output = str(snakemake.output[0])
sample_uid = snakemake.wildcards.sample_uid

# load raw .h5 molecule info from count pipeline

def generate_ir_h5ad(gex_file, vdj_file):
    if "h5ad" in gex_file:
        adata = sc.read_h5ad(gex_file)
    else:
        adata = sc.read_10x_h5(gex_file)
    if "json" in vdj_file:
        adata_bcr = ir.io.read_10x_vdj(vdj_file, include_fields = None)
    if ".tsv" in vdj_file:
        adata_bcr = ir.io.read_airr(vdj_file, include_fields = None)
    ir.pp.merge_airr_chains(adata, adata_bcr)
    adata.obs['sample_uid'] = sample_uid
    return adata

adata_imm = generate_ir_h5ad(gex_file, airr_imm)

adata = generate_ir_h5ad(gex_file, json_TXG)

immcantation = adata_imm.obs

adata.obs = pd.merge(adata.obs, immcantation['IR_VDJ_1_v_identity'], left_index=True, right_index=True, how = 'left')
adata.obs = pd.merge(adata.obs, immcantation['IR_VDJ_1_sequence'], left_index=True, right_index=True, how = 'left')
cut_labels_4 = ['heavily mutated', 'mutated', 'germline']
cut_bins = [0, 0.9, 0.99, 1.000]
adata.obs['mutation_status'] = pd.cut(adata.obs['IR_VDJ_1_v_identity'], bins=cut_bins, labels=cut_labels_4)
# make simple mutation column
cut_labels_2 = ['mutated', 'germline']
cut_bins = [0, 0.99, 1.000]
adata.obs['simple_mutation_status'] = pd.cut(adata.obs['IR_VDJ_1_v_identity'], bins=cut_bins, labels=cut_labels_2)

adata.write_h5ad(merged_output)

print("Done merging VDJ with GEX!")
