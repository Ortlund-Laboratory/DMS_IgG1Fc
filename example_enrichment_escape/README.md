# Example of Enrichment/Escape Snakemake Workflow (FcγR2b)

Fc variant binding to FcγR2b was investigated, and reference, enrichment and escape populations collected for sequencing. These sequences can be processed to give enrichment and escape scores.

## Input Files Required

**SnakeFile**<br>
Gives overall instructions for the `snakemake` workflow.<br>
**config.yaml**<br>
Configuration script controlling variables used by Jupyter notebooks.<br>
[build_variants.ipynb](https://github.com/Ortlund-Laboratory/DMS_IgG1Fc/blob/main/example_enrichment_escape/build_variants.ipynb)<br>
Builds a barcode variant table based on the data from the processed PacBio CCSs.<br>
[R2_to_R1.py](https://github.com/Ortlund-Laboratory/DMS_IgG1Fc/blob/main/example_enrichment_escape/R2_to_R1.py)<br>
Converts barcodes located at the R2 end to the R1 end by taking the reverse complement. This allows the barcodes to be read and parsed correctly by the [illuminabarcodeparser](https://jbloomlab.github.io/dms_variants/dms_variants.illuminabarcodeparser.html#dms_variants.illuminabarcodeparser.IlluminaBarcodeParser) algorithm.<br>
[count_variants.ipynb](https://github.com/Ortlund-Laboratory/DMS_IgG1Fc/blob/main/example_enrichment_escape/count_variants.ipynb)<br>
Counts the number of times a barcode (and by extension a variant) appears in each Illumina barcode sequencing sample.<br>
**counts_to_scores.ipynb**<br>
Groups barcodes and computes enrichment and/or escape scores.<br>
[scripts/run_nb.py](https://github.com/Ortlund-Laboratory/DMS_IgG1Fc/blob/main/example_enrichment_escape/scripts/run_nb.py)<br>
Runs Jupyter notebooks and creates Mardown output.<br>
[data/feature_parse_specs.yaml](https://github.com/Ortlund-Laboratory/DMS_IgG1Fc/blob/main/example_enrichment_escape/data/feature_parse_specs.yaml)<br>
Script for controlling the sequence parsing strategy.<br>
[data/PacBio_amplicons.gb](https://github.com/Ortlund-Laboratory/DMS_IgG1Fc/blob/main/example_enrichment_escape/data/PacBio_amplicons.gb)<br>
GeneBank data file describing sequence features.<br>
[data/barcode_runs.csv](https://github.com/Ortlund-Laboratory/DMS_IgG1Fc/blob/main/example_enrichment_escape/data/barcode_runs.csv)<br>
List of Illumina barcode samples to be analyzed by the snakemake workflow.<br>
[data/processed_ccs.csv](https://github.com/Ortlund-Laboratory/DMS_IgG1Fc/blob/main/example_enrichment_escape/data/processed_ccs.csv)<br>
Processed PacBio CCSs, generated from our [PacBio_Library_Sequencing](https://github.com/Ortlund-Laboratory/DMS_IgG1Fc/tree/main/PacBio_Library_Sequencing) routine. Ensure the library is consistent with those used for the assay.<br>
[data/wildtype_sequence.fasta](https://github.com/Ortlund-Laboratory/DMS_IgG1Fc/blob/main/example_enrichment_escape/data/wildtype_sequence.fasta)<br>
Fc wildtype sequence.<br>


### Sequencing Data

The workflow operates on Illumina barcode sequencing data in fastq.gz format and these files are kept compressed throughout. File location and name should match the listings given in [data/barcode_runs.csv](https://github.com/Ortlund-Laboratory/DMS_IgG1Fc/blob/main/example_enrichment_escape/data/barcode_runs.csv). These files are too large to be contained in GitHub, and so are found, respectively, at:<br>
**p23042-s001_FcLib-1ref_S22_L001_R2_001.fastq.gz** (give link here)<br>
**p23042-s006_JD280top_S27_L001_R2_001.fastq.gz** (give link here)<br>
**p23042-s007_JD280bottom_S28_L001_R2_001.fastq.gz** (give link here)<br>

## Workflow to Count Variants

Use the `snakemake` environment:

`conda activate snakemake`

Run `snakemake` using specified number of cores:

`snakemake -j 6`

## Key Output

[results/counts/barcode_fates.csv](https://github.com/Ortlund-Laboratory/DMS_IgG1Fc/blob/main/example_enrichment_escape/results/counts/barcode_fates.csv)<br>
Tally of barcodes classified and filtered according to quality.<br>
[results/counts/variant_counts.csv](https://github.com/Ortlund-Laboratory/DMS_IgG1Fc/blob/main/example_enrichment_escape/results/counts/variant_counts.csv)<br>
Tally of individual barcode counts for each sample.<br>

## Convert Counts to Enrichment/Escape Scores and Visualize Data

Go to [scores_and_visualization](https://github.com/Ortlund-Laboratory/DMS_IgG1Fc/tree/main/example_enrichment_escape/scores_and_visualization) for count to enrichment/escape score conversions, and data visualization.


