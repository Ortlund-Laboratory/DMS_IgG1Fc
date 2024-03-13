## Input Files Required

[BarcodeMapping_FcgR2b_escape.R](https://github.com/Ortlund-Laboratory/DMS_IgG1Fc/blob/main/example_enrichment_escape/scores_and_visualization/escape/BarcodeMapping_FcgR2b_escape.R)<br>
R script to calculate escape scores from reference and escape counts, generate heatmaps and produce files to map scores onto structures.<br>
[ref_variant_counts.txt](https://github.com/Ortlund-Laboratory/DMS_IgG1Fc/blob/main/example_enrichment_escape/scores_and_visualization/escape/ref_variant_counts.txt)<br>
Barcode counts for the reference sample.<br>
[escape_variant_counts.txt](https://github.com/Ortlund-Laboratory/DMS_IgG1Fc/blob/main/example_enrichment_escape/scores_and_visualization/enrichment/enrich_variant_counts.txt)<br>
Barcode counts for the escape sample.<br>
[Fc_prot.fasta](https://github.com/Ortlund-Laboratory/DMS_IgG1Fc/blob/main/example_enrichment_escape/scores_and_visualization/enrichment/Fc_prot.fasta)
Amino acid sequence for WT Fc. This is required to complete the heatmap.

## Workflow

```
rstudio BarcodeMapping_FcgR2b_escape.R
```

## Key Output

[FcgR2b_Escape_fractions.csv](https://github.com/Ortlund-Laboratory/DMS_IgG1Fc/blob/main/example_enrichment_escape/scores_and_visualization/enrichment/output/FcgR2b_Enrich_fractions.csv)<br>
Log of each mutation and its associated escape score.<br>
[FcgR2b_Escape_Fraction_heatmap01.png](https://github.com/Ortlund-Laboratory/DMS_IgG1Fc/blob/main/example_enrichment_escape/scores_and_visualization/enrichment/output/FcgR2b_Enrich_Fraction_heatmap01.png)<br>
Heatmap of escape scores for all single-point FcγR2b variants.<br>
