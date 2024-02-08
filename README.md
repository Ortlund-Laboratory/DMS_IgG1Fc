# DMS_IgG1Fc
Repo for research related to paper, Engineering IgG1 Fc specificity and affinity for Fc γ receptors using deep mutational scanning titrations
## Paper
Paper can be found [here]need to provide link when available.

Authors: Alasdair D. Keith, Ting Xu, Meredith M. Keen, Tatiana Chernova, Anamika B. Patel, Filipp Frank, Jeffrey V. Ravetch, Eric A. Ortlund and Eric J. Sundberg
## Analysis Workflow
### Build Computing Environment
`conda` is required, and can be obtained _via_ the minimal installer, `miniconda` [here](https://docs.anaconda.com/free/miniconda/).

Installation of mamba (another package manager) is then required since snakemake depends on it. Mamba sometimes will not install with default conda settings. A workaround is to change conda settings to 4.12:

`conda install conda=4.12`
### Input Data

