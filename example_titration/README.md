# Example Titration Snakemake Workflow (FcγR2b)

Fc variant binding to FcγR2b was investigated, with assays run against a range of receptor concentrations, and gated populations then submitted for sequencing. These sequences can be processed to give *K*<sub>D</sub> values.

## Input Files Required

[SnakeFile](https://github.com/Ortlund-Laboratory/DMS_IgG1Fc/blob/main/example_titration/Snakefile)<br>
Gives overall instructions for the `snakemake` workflow.<br>
[config.yaml](https://github.com/Ortlund-Laboratory/DMS_IgG1Fc/blob/main/example_titration/config.yaml)<br>
Configuration script controlling variables used by Jupyter notebooks.<br>
[build_variants.ipynb](https://github.com/Ortlund-Laboratory/DMS_IgG1Fc/blob/main/example_titration/build_variants.ipynb)<br>
Builds a barcode variant table based on the data from the processed PacBio CCSs.<br>
[R2_to_R1.py](https://github.com/Ortlund-Laboratory/DMS_IgG1Fc/blob/main/example_titration/R2_to_R1.py)<br>
Converts barcodes located at the R2 end to the R1 end by taking the reverse complement. This allows the barcodes to be read and parsed correctly by the [illuminabarcodeparser](https://jbloomlab.github.io/dms_variants/dms_variants.illuminabarcodeparser.html#dms_variants.illuminabarcodeparser.IlluminaBarcodeParser) algorithm.<br>
[count_variants.ipynb](https://github.com/Ortlund-Laboratory/DMS_IgG1Fc/blob/main/example_titration/count_variants.ipynb)<br>
Counts the number of times a barcode (and by extension a variant) appears in each Illumina barcode sequencing sample.<br>
[scripts/run_nb.py](https://github.com/Ortlund-Laboratory/DMS_IgG1Fc/blob/main/example_titration/scripts/run_nb.py)<br>
Runs Jupyter notebooks and creates Markdown output.<br>
[data/feature_parse_specs.yaml](https://github.com/Ortlund-Laboratory/DMS_IgG1Fc/blob/main/example_titration/data/feature_parse_specs.yaml)<br>
Script for controlling the sequence parsing strategy.<br>
[data/PacBio_amplicons.gb](https://github.com/Ortlund-Laboratory/DMS_IgG1Fc/blob/main/example_titration/data/PacBio_amplicons.gb)<br>
GeneBank data file describing sequence features.<br>
[data/barcode_runs.csv](https://github.com/Ortlund-Laboratory/DMS_IgG1Fc/blob/main/example_titration/data/barcode_runs.csv)<br>
List of Illumina barcode samples to be analyzed by the snakemake workflow.<br>
[data/processed_ccs.csv](https://github.com/Ortlund-Laboratory/DMS_IgG1Fc/blob/main/example_titration/data/processed_ccs.csv)<br>
Processed PacBio CCSs, generated from our [PacBio_Library_Sequencing](https://github.com/Ortlund-Laboratory/DMS_IgG1Fc/tree/main/PacBio_Library_Sequencing) routine. Ensure the library is consistent with those used for the assay.<br>
[data/wildtype_sequence.fasta](https://github.com/Ortlund-Laboratory/DMS_IgG1Fc/blob/main/example_titration/data/wildtype_sequence.fasta)<br>
Fc wildtype sequence.<br>
[FcgR2b_compute_binding_Kd.Rmd](https://github.com/Ortlund-Laboratory/DMS_IgG1Fc/blob/main/example_titration/FcgR2b_compute_binding_Kd.Rmd)<br>
Converts barcode counts into barcode *K*<sub>D</sub>s. Also, if required, can do weighted averaging of barcodes to give binding curves for selected variants.<br>

### Sequencing Data

The workflow operates on Illumina barcode sequencing data in fastq.gz format and these files are kept compressed throughout. File location and name should match the listings given in **data/barcode_runs.csv**. These files are too large to be contained in GitHub, and so are found, respectively, at:<br>
**p23096-s001_1-1_S212_L002_R2_001.fastq.gz** (give link here)<br>
**p23096-s002_1-2_S213_L002_R2_001.fastq.gz** (give link here)<br>
**p23096-s003_1-3_S214_L002_R2_001.fastq.gz** (give link here)<br>
**p23096-s004_1-4_S215_L002_R2_001.fastq.gz** (give link here)<br>
**p23096-s005_2-1_S216_L002_R2_001.fastq.gz** (give link here)<br>
**p23096-s006_2-2_S217_L002_R2_001.fastq.gz** (give link here)<br>
**p23096-s007_2-3_S218_L002_R2_001.fastq.gz** (give link here)<br>
**p23096-s008_2-4_S219_L002_R2_001.fastq.gz** (give link here)<br>
**p23096-s009_3-1_S220_L002_R2_001.fastq.gz** (give link here)<br>
**p23096-s010_3-2_S221_L002_R2_001.fastq.gz** (give link here)<br>
**p23096-s011_3-3_S222_L002_R2_001.fastq.gz** (give link here)<br>
**p23096-s012_3-4_S223_L002_R2_001.fastq.gz** (give link here)<br>
**p23096-s013_4-1_S224_L002_R2_001.fastq.gz** (give link here)<br>
**p23096-s014_4-2_S225_L002_R2_001.fastq.gz** (give link here)<br>
**p23096-s015_4-3_S226_L002_R2_001.fastq.gz** (give link here)<br>
**p23096-s016_4-4_S227_L002_R2_001.fastq.gz** (give link here)<br>
**p23096-s017_5-1_S228_L002_R2_001.fastq.gz** (give link here)<br>
**p23096-s018_5-2_S229_L002_R2_001.fastq.gz** (give link here)<br>
**p23096-s019_5-3_S230_L002_R2_001.fastq.gz** (give link here)<br>
**p23096-s020_5-4_S231_L002_R2_001.fastq.gz** (give link here)<br>
**p23096-s021_6-1_S232_L002_R2_001.fastq.gz** (give link here)<br>
**p23096-s022_6-2_S233_L002_R2_001.fastq.gz** (give link here)<br>
**p23096-s023_6-3_S234_L002_R2_001.fastq.gz** (give link here)<br>
**p23096-s024_6-4_S235_L002_R2_001.fastq.gz** (give link here)<br>
**p23096-s025_7-1_S236_L002_R2_001.fastq.gz** (give link here)<br>
**p23096-s026_7-2_S237_L002_R2_001.fastq.gz** (give link here)<br>
**p23096-s027_7-3_S238_L002_R2_001.fastq.gz** (give link here)<br>
**p23096-s028_7-4_S239_L002_R2_001.fastq.gz** (give link here)<br>
**p23096-s029_8-1_S240_L002_R2_001.fastq.gz** (give link here)<br>
**p23096-s030_8-2_S241_L002_R2_001.fastq.gz** (give link here)<br>
**p23096-s031_8-3_S242_L002_R2_001.fastq.gz** (give link here)<br>
**p23096-s032_9-1_S243_L002_R2_001.fastq.gz** (give link here)<br>
**p23096-s033_9-2_S244_L002_R2_001.fastq.gz** (give link here)<br>
**p23096-s034_9-3_S245_L002_R2_001.fastq.gz** (give link here)<br>
**p23096-s035_10-1_S246_L002_R2_001.fastq.gz** (give link here)<br>
**p23096-s036_10-2_S247_L002_R2_001.fastq.gz** (give link here)<br>
**p23096-s037_10-3_S248_L002_R2_001.fastq.gz** (give link here)<br>
**p23096-s038_11-1_S249_L002_R2_001.fastq.gz** (give link here)<br>
**p23096-s039_11-2_S250_L002_R2_001.fastq.gz** (give link here)<br>
**p23096-s040_11-3_S251_L002_R2_001.fastq.gz** (give link here)<br>

**NOTE**: Not all bins produced a large enough population for sequencing. However, the analysis software requires all bins to be present, and dms_variants cannot read empty files. Therefore, for each empty bin, a filler file prefixed *sparefile_*, is added. These files are all formatted as follows:

```
@Sparefile 2:N:0
GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
+
FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
```
This will ensure the code runs smoothly but will not affect analysis since there is no sequence in the library which contains so many sequential Gs.

In this example, the filler files required were:<br>
[sparefile_conc0_bin4.fastq.gz](https://github.com/Ortlund-Laboratory/DMS_IgG1Fc/blob/main/example_titration/sparefile_conc0_bin4.fastq.gz)<br>
[sparefile_conc1_bin4.fastq.gz](https://github.com/Ortlund-Laboratory/DMS_IgG1Fc/blob/main/example_titration/sparefile_conc1_bin4.fastq.gz)<br>
[sparefile_conc2_bin4.fastq.gz](https://github.com/Ortlund-Laboratory/DMS_IgG1Fc/blob/main/example_titration/sparefile_conc2_bin4.fastq.gz)<br>
[sparefile_conc3_bin4.fastq.gz](https://github.com/Ortlund-Laboratory/DMS_IgG1Fc/blob/main/example_titration/sparefile_conc3_bin4.fastq.gz)<br>

## Snakemake Workflow

Use the `snakemake` environment:

`conda activate snakemake`

Run `snakemake` using specified number of cores:

`snakemake -j 6`

## Snakemake Key Output

[results/counts/barcode_fates.csv](https://github.com/Ortlund-Laboratory/DMS_IgG1Fc/blob/main/example_titration/results/counts/barcode_fates.csv)<br>
Tally of barcodes classified and filtered according to quality.<br>
**results/counts/variant_counts.csv** (give link here)<br>
Tally of individual barcode counts for each sample.<br>

## *K*<sub>D</sub> Generation Workflow

In the same directory, run:

```
rstudio FcgR2b_compute_binding_Kd.Rmd
```

## *K*<sub>D</sub> Generation Key Output

[binding_Kds.csv](https://github.com/Ortlund-Laboratory/DMS_IgG1Fc/blob/main/example_titration/results/binding_Kds/binding_Kds.csv)<br>
Log of calculated *K*<sub>D</sub> values for each barcode.<br>

## Grouping Barcodes Workflow



## Data Visualization Workflow


