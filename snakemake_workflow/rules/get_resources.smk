rule get_transcriptome_ref:
    output:
        transcriptome="transcriptome".format(workflow.basedir),
        outdir=directory("{}/db/transcriptome/".format(workflow.basedir)),
    params:
        ref="{}/db/10X/{}/".format(workflow.basedir, config["vdj_ref_prefix"]),
        outdir="{}/db/10X/".format(workflow.basedir),
        species="Homo sapiens",
        link=config["10X_transcriptome"],
        name=config["10X_transcriptome"].split("/")[-1]
    conda:
        "../envs/cellranger.yaml"
    shell:
        "mkdir -p {output.outdir} "
        "&& cd {output.outdir} "
        "&& wget {params.link} "
        "&& tar -xvf {params.name}"
