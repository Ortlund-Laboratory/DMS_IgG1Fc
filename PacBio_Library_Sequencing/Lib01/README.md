# PacBio Library Sequencing - Workflow

## Input Files Required

### config.yaml
Configuration script controlling variables used by Jupyter notebook.
### data/feature_parse_specs.yaml
Script for controlling the sequence parsing strategy.
### data/PacBio_amplicons.gb
GeneBank data file describing sequence features.
### data/PacBio_runs.csv
List of sequence (fastq) files to be analyzed.
### results/ccs/XXX.fastq
Input circular consensus sequences (CCSs) data file.
### process_ccs.ipynb
Jupyter notebook for extracting barcodes from CCSs and matching to variants.

## Setup

Ensure the fastq file(s) listed in PacBio_runs.csv are present in results/ccs BEFORE executing the workflow.

## Workflow

Use the 
