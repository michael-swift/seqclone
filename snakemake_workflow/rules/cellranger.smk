import sys
import os

cellbender_prior_cells = 6000

def expected_cells_adj(wildcards, attempt):
    return cellbender_prior_cells / attempt

#### CellRanger VDJ #####
multi_config = "multi_config"

rule cellranger_multi:
    input:
        workflow.basedir + "/config/multi_config_{sample_uid}.csv",
    output:directory("{base}/per_sample/txg/{sample_uid}"),
        h5="{base}/per_sample/txg/{sample_uid}/outs/multi/count/raw_feature_bc_matrix.h5",
        airr="{base}/per_sample/txg/{sample_uid}/outs/per_sample_outs/{sample_uid}/vdj_b/airr_rearrangement.tsv",
        json="{base}/per_sample/txg/{sample_uid}/outs/multi/vdj_b/all_contig_annotations.json",
        contigs="{base}/per_sample/txg/{sample_uid}/outs/multi/vdj_b/all_contig.fasta",
        anno="{base}/per_sample/txg/{sample_uid}/outs/multi/vdj_b/all_contig_annotations.csv",
    log:
        "{base}/logs/{sample_uid}/cellranger.log",
    params:
        name="count_cellranger",
        base=config["base"],
        cell_ranger=config["cell_ranger"],
    resources:
        mem_mb=210000,
    threads: 20
    shell:
        "mkdir -p {params.base}/txg &&"
        " cd {params.base}/txg &&"
        " rm -rf {wildcards.sample_uid} && {params.cell_ranger}/cellranger multi --id={wildcards.sample_uid} --csv={input}"

rule run_cellbender:
    input:
        ancient("{base}/per_sample/txg/{sample_uid}/outs/multi/count/raw_feature_bc_matrix.h5")
    output:
        "{base}/per_sample/cellbender/{sample_uid}/background_removed.h5"
    conda:
        "cellbender2"
    resources:
        expected_cells=expected_cells_adj
    log:
        "{base}/logs/{sample_uid}/cellbender.log"
    shell:
        "cellbender remove-background"
        " --input {input}"
        " --output {output}"
        " --expected-cells {resources.expected_cells}"
        " --fpr 0.01"
        " --cuda"

rule create_immcantation_dirs:
    input:
        contigs = ancient("{base}/per_sample/txg/{sample_uid}/outs/multi/vdj_b/all_contig.fasta"),
        anno = ancient("{base}/per_sample/txg/{sample_uid}/outs/multi/vdj_b/all_contig_annotations.csv")
    output:
        anno="{base}/per_sample/immcantation/{sample_uid}/all_contig_annotations.csv",
        contigs="{base}/per_sample/immcantation/{sample_uid}/all_contig.fasta",
    resources:
        mem_mb=8000,
        disk_mb=4000,
    params:
        immcantation_dir="{base}/per_sample/immcantation/{sample_uid}",
    shell:
        "mkdir -p {params.immcantation_dir} "
        " && cp {input.contigs} {params.immcantation_dir}"
        " && cp {input.anno} {params.immcantation_dir}"

rule immcantation_changeo10x:
    input:
        contigs=rules.create_immcantation_dirs.output.contigs,
        anno=rules.create_immcantation_dirs.output.anno,
    output:
        immcantation="{base}/per_sample/immcantation/{sample_uid}/{sample_uid}_heavy_germ-pass.tsv",
    log:
        "{base}/logs/{sample_uid}/immcantation.log",
    params:
        name="immcantation",
        base=config["base"],
        sif="{}/{}".format(config["resources"], "immcantation_suite-4.3.0.sif"),
        immcantation_dir="{base}/per_sample/immcantation/{sample_uid}",
    shell:
        "module load system && "
        "cd {params.immcantation_dir} && "
        "singularity exec -B {params.immcantation_dir}:/data {params.sif} "
        "changeo-10x -s {input.contigs} -a {input.anno} -f airr -x 0.10 -n {wildcards.sample_uid} -o {params.immcantation_dir}"

rule scirpy_merge:
    input:
        h5="{base}/per_sample/cellbender/{sample_uid}/background_removed.h5",
        airr_imm="{base}/per_sample/immcantation/{sample_uid}/{sample_uid}_heavy_germ-pass.tsv",
        json_txg=ancient("{base}/per_sample/txg/{sample_uid}/outs/multi/vdj_b/all_contig_annotations.json"),
    output:
        "{base}/per_sample/scirpy/{sample_uid}/ir_object.h5ad"
    resources:
        mem_mb=16000,
        partition="quake,owners",
        time="0-1",
    log:
        "{base}/logs/per_sample/{sample_uid}/scirpy.log",
    params:
        name="scirpy",
        base=config["base"]
    conda:"/home/groups/quake/mswift/envs/scanpy_latest"
#          "scanpy_latest"
    script:
        config["workflow_dir"] + "scripts/scirpy_merge.py"

rule aggregate_gex:
    """ Aggregates all the genes expression data from txg with very lenient filtering """
    input:
        expand(
            "{base}/per_sample/scirpy/{sample_uid}/ir_object.h5ad",
            base=config["base"],
            sample_uid=sample_uids),
    output:
        temp("{base}/aggregated/merged/aggr_gex.h5ad"),
    log:
        "{base}/logs/merged_aggregation.log",
    resources:
        mem_mb=32000,
        partition="quake,owners",
        time="0-1",
    conda:
        config["workflow_dir"] + "envs/scirpy.yaml"
    params:
        min_genes="200",
        min_counts="500",
    script:
        config["workflow_dir"] + "scripts/aggregate_with_scanpy.py"

rule preprocess_scirpy:
    input:
        "{base}/aggregated/merged/aggr_gex.h5ad",
    output:
        "{base}/processed/merged/scirpy_processed.h5ad",
    log:
        "{base}/logs/merged_preprocess.log",
    resources:
        mem_mb=64000,
        partition="quake",
        time="0-1",
    conda:
        config["workflow_dir"] + "envs/scirpy.yaml"
    params:
        #filter leniently just for memory reasons
        min_genes="200",
        min_counts="500",
    script:
        config["workflow_dir"] + "scripts/assign_clones_scirpy.py"
