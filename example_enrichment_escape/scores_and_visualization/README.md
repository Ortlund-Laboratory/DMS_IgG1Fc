# File Prep

[variant_counts.csv](https://github.com/Ortlund-Laboratory/DMS_IgG1Fc/blob/main/example_enrichment_escape/results/counts/variant_counts.csv) from results/counts needs to be manipulated.

In our example, by referring to [barcode_runs.csv](https://github.com/Ortlund-Laboratory/DMS_IgG1Fc/blob/main/example_enrichment_escape/data/barcode_runs.csv), we can see that three experiments were run: a reference sample (exp01-none-0-reference), an  enrichment sample with FcγR2b (exp02-JD280top-2000-escape) and an escape sample with FcγR2b (exp03-JD280bottom-2000-escape). 

**NOTE**: Though the above enrichment sample has escape in its title, this was due to naming conventions in some of the Jupyter scripts, which requires samples to be classified as either 'reference' or 'escape'. In our naming system, samples which have 'top' in their titles denote enrichment populations, and samples with 'bottom' in their titles denote escape populations.

For compatibility with our R scripts, we require the information in [variant_counts.csv](https://github.com/Ortlund-Laboratory/DMS_IgG1Fc/blob/main/example_enrichment_escape/results/counts/variant_counts.csv) to be split into three-column (barcode,mutation,count) files.

First, remove all variants where the number of amino acid mutations is 0 (i.e. is still wildtype) or greater than 1 (multiple point mutations).

```
awk -F',' '!($10!=1)' variant_counts.csv > tmp.csv && mv tmp.csv variant_counts.csv
```

```
sed -n '/exp01-none-0-reference/p' variant_counts.csv > ref_variant_counts.txt
sed -n '/exp02-JD280top-2000-escape/p' variant_counts.csv > enrich_variant_counts.txt
sed -n '/exp03-JD280bottom-2000-escape/p' variant_counts.csv > escape_variant_counts.txt
```

