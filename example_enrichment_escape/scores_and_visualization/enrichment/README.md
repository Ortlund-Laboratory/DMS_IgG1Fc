## Input Files Required

[BarcodeMapping_FcgR2b_enrich.R](https://github.com/Ortlund-Laboratory/DMS_IgG1Fc/blob/main/example_enrichment_escape/scores_and_visualization/enrichment/BarcodeMapping_FcgR2b_enrich.R)<br>
R script to calculate enrichment scores from reference and enrichment counts, generate heatmaps and produce files to map scores onto structures.<br>
[ref_variant_counts.txt](https://github.com/Ortlund-Laboratory/DMS_IgG1Fc/blob/main/example_enrichment_escape/scores_and_visualization/enrichment/ref_variant_counts.txt)<br>
Barcode counts for the reference sample.<br>
[enrich_variant_counts.txt](https://github.com/Ortlund-Laboratory/DMS_IgG1Fc/blob/main/example_enrichment_escape/scores_and_visualization/enrichment/enrich_variant_counts.txt)<br>
Barcode counts for the enrichment sample.<br>
[Fc_prot.fasta](https://github.com/Ortlund-Laboratory/DMS_IgG1Fc/blob/main/example_enrichment_escape/scores_and_visualization/enrichment/Fc_prot.fasta)
Amino acid sequence for WT Fc. This is required to complete the heatmap.

## Workflow

```
rstudio BarcodeMapping_FcgR2b_escape.R
```
