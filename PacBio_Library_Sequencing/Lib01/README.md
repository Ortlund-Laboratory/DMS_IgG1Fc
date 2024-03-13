# PacBio Library Sequencing - Workflow

## Input Files Required

[config.yaml](https://github.com/Ortlund-Laboratory/DMS_IgG1Fc/blob/main/PacBio_Library_Sequencing/Lib01/config.yaml)<br>
Configuration script controlling variables used by Jupyter notebook.<br>
[data/feature_parse_specs.yaml](https://github.com/Ortlund-Laboratory/DMS_IgG1Fc/blob/main/PacBio_Library_Sequencing/Lib01/data/feature_parse_specs.yaml)<br>
Script for controlling the sequence parsing strategy.<br>
[data/PacBio_amplicons.gb](https://github.com/Ortlund-Laboratory/DMS_IgG1Fc/blob/main/PacBio_Library_Sequencing/Lib01/data/PacBio_amplicons.gb)<br>
GeneBank data file describing sequence features.<br>
**data/PacBio_runs.csv**<br>
List of sequence (fastq) files to be analyzed.<br>
**results/ccs/XXX.fastq**<br>
Input circular consensus sequences (CCSs) data file.<br>
**process_ccs.ipynb**<br>
Jupyter notebook for extracting barcodes from CCSs and matching to variants.<br>

## Setup

Ensure the fastq file(s) listed in PacBio_runs.csv are present in results/ccs BEFORE executing the workflow.

## Workflow

Use the `snakemake` environment:

`conda activate snakemake`

Run `jupyter`:

`jupyter notebook process_ccs.ipynb`

**NOTE:** Some cells of the `jupyter` notebook may fail to execute, in the absence of ccs summary files. These are not required, and so the user can skip to the next cell.

## Key Output

**results/process_ccs/processed_ccs.csv**<br>
List of parsed barcodes and their associated mutations.

