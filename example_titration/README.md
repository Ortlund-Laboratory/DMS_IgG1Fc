# Example of Titration Snakemake Workflow (FcγR2b)

Fc variant binding to FcγR2b was investigated, with assays run against a range of receptor concentrations, and gated populations then submitted for sequencing. These sequences can be processed to give *K*<sub>D</sub> values.

## Input Files Required

**SnakeFile**<br>
Gives overall instructions for the `snakemake` workflow.<br>
**config.yaml**<br>
Configuration script controlling variables used by Jupyter notebooks.<br>
**build_variants.ipynb**<br>
Builds a barcode variant table based on the data from the processed PacBio CCSs.<br>
**R2_to_R1.py**<br>
Converts barcodes located at the R2 end to the R1 end by taking the reverse complement. This allows the barcodes to be read and parsed correctly by the [illuminabarcodeparser](https://jbloomlab.github.io/dms_variants/dms_variants.illuminabarcodeparser.html#dms_variants.illuminabarcodeparser.IlluminaBarcodeParser) algorithm.<br>
**count_variants.ipynb**<br>
Counts the number of times a barcode (and by extension a variant) appears in each Illumina barcode sequencing sample.<br>
**scripts/run_nb.py**<br>
Runs Jupyter notebooks and creates Mardown output.<br>
**data/feature_parse_specs.yaml**<br>
Script for controlling the sequence parsing strategy.<br>
**data/PacBio_amplicons.gb**<br>
GeneBank data file describing sequence features.<br>
**data/barcode_runs.csv**<br>
List of Illumina barcode samples to be analyzed by the snakemake workflow.<br>
**data/processed_ccs.csv**<br>
Processed PacBio CCSs, generated from our [PacBio_Library_Sequencing](https://github.com/Ortlund-Laboratory/DMS_IgG1Fc/tree/main/PacBio_Library_Sequencing) routine. Ensure the library is consistent with those used for the assay.<br>
**data/wildtype_sequence.fasta**<br>
Fc wildtype sequence.<br>
