import os
import pandas as pd
import glob
from collections import defaultdict

shell.prefix("set +euo pipefail;")


configfile: "config/config.yaml"


sample_uids = config["10X_lanes"]
base = config["base"]
# should be loaded as config file

os.makedirs(base, exist_ok=True)


rule all:
    input:
        expand("{base}/processed/merged/scirpy_processed.h5ad", base = config["base"])
    params:
        name="all",
        partition="quake",
    threads: 1


include: "rules/get_resources.smk"
include: "rules/cellranger.smk"

localrules: all
