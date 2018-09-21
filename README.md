# DM.R
### A workflow for the statistical analysis and visualization of differentially methylated regions (DMRs) from a CpG count matrix

## Installation

No manual installation of R packages is required, the required packages and updates will occur automaitcally upon running the script.

## The Design Matrix and Covariates

This script requires a basic design matrix to identify the groups and covariates, which should be named `sample_info.csv` and contain header columns to idenfity the factor. It is important to have an experimental sample, rather than a control sample, as the first sample in the design matrix in order to obtain results for experimental vs. control rather than control vs. experiemntal.

Within the script, covariates can be selected for adjustment. There are two ways to adjust for covariates:
1. Use `adjustCovariate` for covariates that are continuous or contain two or more groups. More than one covariate can be adjusted for.
2. Use `matchCovariate` for balancing permutations, which is ideal for two group covariates such as sex. Only one covariate can be balanced.

## Input

This workflow requires the following variables:
1. genome: Select either: hg38, mm10, rn6, or rheMac8
2. coverage: The coverage cutoff for all samples, 1x is recommended
3. testCovariate: The covariate to test for signficant differences between experimental and control, i.e.: Diagnosis
4. adjustCovariate: The covariates to directly adjust, i.e.: Age
5. matchCovariate: The covariate to balance permutations, i.e: Sex
6. cores: The number of cores to use, 2 are recommended

## Output

This workflow provides the following files:
1. CpG methylation and coverage value distrubution plots
2. DMRs and background regions
3. Individual smoothed methylation values for DMRs, background regions, 20 Kb windows, and CpG island windows
4. Heatmap of DMRs
5. PCA plots of 20 Kb windows and CpG island windows
6. Gene ontologies and pathways (enrichr for all, GREAT for hg38)
7. Gene region and CpG annotations and plots (hg38, mm10, or rn6)
8. Manhattan and Q-Qplots 
9. Blocks of methylation and background blocks

## Conceptual Questions

This workflow is primarily based on the [dmrseq](https://www.bioconductor.org/packages/release/bioc/html/dmrseq.html) and [bsseq](https://www.bioconductor.org/packages/release/bioc/html/bsseq.html) bioconductor packages.

While developing this script, I was fortunate enough to have [Keegan Korthauer](https://github.com/kdkorthauer), the creator of dmrseq, provide invaluable insight:

> Keegan re beta coefficent: You’re exactly right that it represents the average effect size over the region, but if you’d like to take it a step further and connect it to the difference seen in the plot, you can divide the beta coefficient by pi (yep, 3.14159…) to put it on the scale of a proportion difference. This is because the beta coefficient is on the scale of the arcsine transformed differences. So beta/pi will be similar to (and correlated with) the simple mean proportion difference across the region, but the beta/pi quantity from the model is adjusted for things like coverage and correlated errors. 

> Keegan re individual values: The per sample smoothing lines in the plots (1) are very different than the smoothed methylation differences dmrseq computes,and (2) are *purely* for visualization purposes. They simply smooth the methylation values with loess, and do not use the model in any way. If you really need smoothed sample-specific methylation values, I’d suggest obtaining them with the bsseq package.

## Acknowledgements

The author would like to thank [Keegan Korthauer](https://github.com/kdkorthauer) for helpful conceptual advice in establishing and optimizing this workflow. The author would also like to thank [Nikhil Joshi](https://github.com/najoshi) from the [UC Davis Bioinformatics Core](https://github.com/ucdavis-bioinformatics) for troubleshooting of a [resource issue](https://github.com/kdkorthauer/dmrseq/commit/38dea275bb53fcff3a0df93895af759b15c90e3e).
