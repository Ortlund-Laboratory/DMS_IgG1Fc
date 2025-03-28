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

```
sort -t',' -k 4 -g FcgR2b-FcgR2a-131H_Ka_fractions.csv > sorted_FcgR2b-FcgR2a-131H_Ka_fractions.csv
```

## Key Output

[FcgR2b-FcgR2a-131H_Ka_fractions.csv](https://github.com/Ortlund-Laboratory/DMS_IgG1Fc/blob/main/example_differential_analysis/FcgR2b-FcgR2a-131H_Ka_fractions.csv)<br>
Differential FcgR2b-FcgR2a-131H log*K*<sub>a</sub> values for each Fc mutation.<br>
[sorted_FcgR2b-FcgR2a-131H_Ka_fractions.csv](https://github.com/Ortlund-Laboratory/DMS_IgG1Fc/blob/main/example_differential_analysis/sorted_FcgR2b-FcgR2a-131H_Ka_fractions.csv)<br>
Differential FcgR2b-FcgR2a-131H log*K*<sub>a</sub> values for each Fc mutation, organized from most FcgR2a-131H favoring to most FcgR2b favoring.<br>
[FcgR2b-FcgR2a-131H_average_Ka.csv](https://github.com/Ortlund-Laboratory/DMS_IgG1Fc/blob/main/example_differential_analysis/FcgR2b-FcgR2a-131H_average_Ka.csv)<br>
Differential FcgR2b-FcgR2a-131H log*K*<sub>a</sub> values averaged over each Fc site.<br>
[FcgR2b-FcgR2a-131H_Normalized_LogKa_Heatmap01.png](https://github.com/Ortlund-Laboratory/DMS_IgG1Fc/blob/main/example_differential_analysis/FcgR2b-FcgR2a-131H_Normalized_LogKa_Heatmap01.png)<br>
Heatmap of differential FcgR2b-FcgR2a-131H log*K*<sub>a</sub> values for each Fc mutation. Green points denote mutations which favor FcgR2b over FcgR2a-131H compared to WT, and pink points denote mutations which favor FcgR2a-131H.<br>
[FcgR2b-FcgR2a-131H_color_average_Ka.png](https://github.com/Ortlund-Laboratory/DMS_IgG1Fc/blob/main/example_differential_analysis/FcgR2b-FcgR2a-131H_color_average_Ka.png)<br>
Bar chart of differential FcgR2b-FcgR2a-131H log*K*<sub>a</sub> values, averaged over each Fc site. Color-coding is same as for the heatmap.<br>
