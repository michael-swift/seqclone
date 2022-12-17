import pandas as pd
import scanpy as sc
import scirpy as ir

# I/O
gex_file = snakemake.input['h5']
airr_imm = snakemake.input['airr_imm']
airr_10X = snakemake.input['airr_10X']
output_imm = str(snakemake.output[1])
output_10X = str(snakemake.output[0])
print(output_imm)
print(output_10X)
# load raw .h5 molecule info from count pipeline
def generate_ir_h5ad(gex_file, airr_file, output_name):
    adata = sc.read_10x_h5(gex_file), cache = False)
    adata_bcr = ir.io.read_airr(airr_file)
    ir.pp.merge_airr_chains(adata, adata_bcr)
    adata.write_h5ad(output_name)
    print("wrote {}".format(output_name))
    return

generate_ir_h5ad(gex_file, airr_imm, output_imm)

generate_ir_h5ad(gex_file, airr_10X, output_10X)

print("Done merging VDJ with GEX!")
