# seqclone

seqclone is an algorithm that identifies gene expression states inherited by clones during differentiation. These genes could become targets for further functional studies and have been useful for quantifying the extent of clonal effects in differentiation. I used it for the analyses in the paper ["Lineage tracing reveals fate bias and transcriptional memory in human B cells"](https://www.life-science-alliance.org/content/lsa/6/3/e202201792.full.pdf), for which you can find the processed data [here](link-to-processed-data).

## Approach

It implements the approach described in ["Multiplexed Division Tracking Dyes for Proliferation-Based Clonal Lineage Tracing"](https://doi.org/10.4049/jimmunol.1800481) but does so in Python, on a transcriptomic scale, using single-cell RNA-seq data (although any other data type in AnnData format with clonal labels should work).

## Usage
The script `run_seqclone.py` is used to identify genes with clonally inherited expression states from single-cell RNA-seq data. Here's a breakdown of its usage:

1. The script takes two command-line arguments:
   - `-config`: Path to the configuration file (default: 'config')
   - `-outdir`: Output directory for the results (default: current working directory)

2. The configuration file is an `.ini` file that contains parameters for filtering and input/output settings. The script reads these parameters using the `configparser` module.

3. The script loads the single-cell RNA-seq data from an `h5ad` file specified in the configuration file using the `scanpy` library.
4. 
To use the script, you need to:
1. Prepare your single-cell RNA-seq data in `h5ad` format with a column in `adata.obs` named "clone_id" indicating the clonal membership of each cell.
2. Create a configuration file specifying the input data, filtering parameters, and output settings.
3. Run the script with the appropriate command-line arguments, specifying the path to the configuration file and the desired output directory.

The script will perform permutation testing to identify genes with significant clonal expression patterns and generate output files with the results.

Note: The script assumes the presence of a "clone_id" column in `adata.obs` and may require modifications to work with different data formats or column names.

## Disclaimer

Please note that this code is not actively maintained. However, I am available to assist with any questions or issues you may encounter while using seqclone. Feel free to contact me for support.


