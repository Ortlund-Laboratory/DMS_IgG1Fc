# Example of Enrichment/Escape Snakemake Workflow (FcγR2b)

Fc variant binding to FcγR2b was investigated, and reference, enrichment and escape populations collected for sequencing. These sequences can be processed to give enrichment and escape scores.

## Input Files Required

**Snakefile**<br>
Gives overall instructions for the `snakemake` workflow.<br>
[build_variants.ipynb](https://github.com/Ortlund-Laboratory/DMS_IgG1Fc/blob/main/example_enrichment_escape/build_variants.ipynb)<br>
Builds a barcode variant table based on the data from the processed PacBio CCSs.<br>
[R2_to_R1.py](https://github.com/Ortlund-Laboratory/DMS_IgG1Fc/blob/main/example_enrichment_escape/R2_to_R1.py)<br>
Converts barcodes located at the R2 end to the R1 end by taking the reverse complement. This allows the barcodes to be read and parsed correctly by the [illuminabarcodeparser](https://jbloomlab.github.io/dms_variants/dms_variants.illuminabarcodeparser.html#dms_variants.illuminabarcodeparser.IlluminaBarcodeParser) algorithm.<br>
[count_variants.ipynb](https://github.com/Ortlund-Laboratory/DMS_IgG1Fc/blob/main/example_enrichment_escape/count_variants.ipynb)<br>
Counts the number of times a barcode (and by extension a variant) appears in each Illumina barcode sequencing sample.<br>
[scripts/run_nb.py](https://github.com/Ortlund-Laboratory/DMS_IgG1Fc/blob/main/example_enrichment_escape/scripts/run_nb.py)<br>
Runs Jupyter notebooks and creates Mardown output.<br>
[data/feature_parse_specs.yaml](https://github.com/Ortlund-Laboratory/DMS_IgG1Fc/blob/main/example_enrichment_escape/data/feature_parse_specs.yaml)<br>
Script for controlling the sequence parsing strategy.<br>
[data/PacBio_amplicons.gb](https://github.com/Ortlund-Laboratory/DMS_IgG1Fc/blob/main/example_enrichment_escape/data/PacBio_amplicons.gb)<br>
GeneBank data file describing sequence features.<br>


### Sequencing Data


## Workflow

Use the `snakemake` environment:

`conda activate snakemake`

Run `snakemake` using specified number of cores:

`snakemake -j 6`
