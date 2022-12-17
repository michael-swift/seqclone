#!/bin/bash
#name shell variables for calling in snakemake and sbatch
DATETIME=$(date "+%Y_%m_%d_%H_%M_%S")
RESTART=0

SNAKEFILE=Snakefile.smk
EMAIL=mswift2@stanford.edu

CLUSTER_CONFIG=config/slurm_config.yaml
#Snakemake config
NJOBS=200
WAIT=120

source /home/groups/quake/mswift/miniconda3/etc/profile.d/conda.sh
conda activate snakemake
mkdir -p snakemake_logs/slurm_logs/

#log file for process that calls snakemake
SBATCH_LOGFILE=snakemake_logs/cluster.$DATETIME.log
SBATCH_LOGFILE_ERR=snakemake_logs/cluster.$DATETIME.log.err


# for a specific cell / output
TARGET=''


if [ $# -eq 0 ]
  then
    # Dry run snakemake
    snakemake -s $SNAKEFILE $TARGET --use-conda --keep-target-files --rerun-incomplete -n -r --quiet --keep-going
elif [ $1 = "unlock" ]
    then
        snakemake -s $SNAKEFILE $TARGET -F --rerun-incomplete --unlock --cores 1
    
elif [ $1 = "dryrun" ]
    # Dry run snakemake and print shell cmds 
    then
        snakemake \
            -s $SNAKEFILE $TARGET \
            --use-conda \
            --keep-target-files \
            --rerun-incomplete \
            -n \
            -r \
            -k \
            --printshellcmds

elif [ $1 = "profiledry" ]
    then
  # Run snakemake
    echo 'running snakemake'
    snakemake --profile profile/ -n

elif [ $1 = "profile" ]
    then
  # Run snakemake
    echo 'running snakemake'
    snakemake --profile profile/

elif [ $1 = "sbatch" ]
    # Run snakemake as an SBATCH job, good for long workflows
    then
    sbatch \
        --ntasks=1 \
        --cpus-per-task=1 \
        --mem=8000 \
        --mail-user=$EMAIL \
        --time 4-0 \
        -p quake \
        -o $SBATCH_LOGFILE \
        -e $SBATCH_LOGFILE_ERR \
        run_snake.sh profile

else
    echo "wrong option"
fi
