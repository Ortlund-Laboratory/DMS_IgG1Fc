# Example Differential Analysis Workflow

Using the processed titration data, it is possible to determine how different receptors compare to one another with regard to binding to the Fc mutants. The example workflow here compares the FcγR2b and FcγR2a-131H receptors.

## Input Files Required

[FcgR2b_mutation_to_Ka.csv](https://github.com/Ortlund-Laboratory/DMS_IgG1Fc/blob/main/example_differential_analysis/FcgR2b_mutation_to_Ka.csv)<br>
This is the output from our FcgR2b titration analysis, renamed from **mutation_to_Ka.csv** to differentiate it from FcgR2a-131H.<br>
[FcgR2a-131H_mutation_to_Ka.csv](https://github.com/Ortlund-Laboratory/DMS_IgG1Fc/blob/main/example_differential_analysis/FcgR2a-131H_mutation_to_Ka.csv)<br>
This is the output from our FcgR2a-131H titration analysis, renamed from **mutation_to_Ka.csv** to differentiate it from FcgR2b.<br>
[FcgR2b-FcgR2a-131H_differential_analysis.R](https://github.com/Ortlund-Laboratory/DMS_IgG1Fc/blob/main/example_differential_analysis/FcgR2b-FcgR2a-131H_differential_analysis.R)<br>
Script for generating the difference in log*K*<sub>a</sub> values for each Fc mutation with respect to the two receptors being compared. Ensure that the log*K*<sub>a</sub> values for WT Fc binding to each receptor is correctly inputted as this is important for normalization.<br>
[Fc_prot.fasta](https://github.com/Ortlund-Laboratory/DMS_IgG1Fc/blob/main/example_differential_analysis/Fc_prot.fasta)<br>
Amino acid sequence for WT Fc. This is required to complete the heatmap.<br>

## Workflow

```
rstudio FcgR2b-FcgR2a-131H_differential_analysis.R
```
Then, if you wish to organize your results from most FcgR2a-131H favoring to FcgR2b favoring:
sort -t',' -k 4 -g FcgR2b-FcgR2a-131H_Ka_fractions.csv > sorted_FcgR2b-FcgR2a-131H_Ka_fractions.csv
