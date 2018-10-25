# DMRichR
### A workflow for the statistical analysis and visualization of differentially methylated regions (DMRs) from a CpG count matrix

## Installation

No manual installation of R packages is required, the required packages and updates will occur automatically upon running the executable script located in the exec folder.

## The Design Matrix and Covariates

This script requires a basic design matrix to identify the groups and covariates, which should be named `sample_info.csv` and contain header columns to identify the factor. It is important to have the label for the experimental samples start with a letter in the alphabet that comes after the one used for control samples in order to obtain results for experimental vs. control rather than control vs. experimental. Within the script, covariates can be selected for adjustment. There are two different ways to adjust for covariates: directly adjust values or balance permutations.


| Name          | Diagnosis      | Age           |  Sex          |
| ------------- | -------------- | ------------- | ------------- |
| SRR3537014    | Idiopathic_ASD | 14            | M             |
| SRR3536981    | Control        | 42            | F             |


## Input

This workflow requires the following variables:
1. `-g --genome` Select either: hg38, mm10, rn6, or rheMac8
2. `-x --coverage` The coverage cutoff for all samples, 1x is recommended
3. `-t --testCovariate` The covariate to test for significant differences between experimental and control, i.e.: Diagnosis
4. `-a --adjustCovariate` Adjust covariates that are continuous or contain two or more groups. More than one covariate can be adjusted for., i.e.: "Age" or c("Age", "PMI")
5. `-m --matchCovariate` Covariate to balance permutations, which is ideal for two group covariates. Only one covariate can be balanced. i.e: Sex
6. `-c --cores` The number of cores to use, 2 are recommended

## Output

This workflow provides the following files:
1. CpG methylation and coverage value distribution plots
2. DMRs and background regions
3. Individual smoothed methylation values for DMRs, background regions, and windows/bins
4. Smoothed global and chromosomal methylation values and statistics
5. Heatmap of DMRs
6. PCA plots of 20 Kb windows (all genomes and hg38, mm10, and rn6 for CpG island windows)
7. Gene ontologies and pathways (enrichr for all genomes, GREAT for hg38 and mm10)
8. Gene region and CpG annotations and plots (hg38, mm10, or rn6)
9. Manhattan and Q-Qplots 
10. Blocks of methylation and background blocks

## Acknowledgements

This workflow is primarily based on the [dmrseq](https://www.bioconductor.org/packages/release/bioc/html/dmrseq.html) and [bsseq](https://www.bioconductor.org/packages/release/bioc/html/bsseq.html) bioconductor packages. I would like to thank [Keegan Korthauer](https://github.com/kdkorthauer), the creator of dmrseq, for helpful conceptual advice in establishing and optimizing this workflow. I would like to thank [Matt Settles](https://github.com/msettles) from the [UC Davis Bioinformatics Core](https://github.com/ucdavis-bioinformatics) for advice on creating an R package and use of the tidyverse. I would like to thank Rochelle Coulson for a script that was developed into the PCA function. I would also like to thank Blythe Durbin-Johnson and Annie Vogel Ciernia for statistical consulting that enabled the global and chromosomal methylation statistics. Finally, I would like to thank [Nikhil Joshi](https://github.com/najoshi) from the [UC Davis Bioinformatics Core](https://github.com/ucdavis-bioinformatics) for troubleshooting of a [resource issue](https://github.com/kdkorthauer/dmrseq/commit/38dea275bb53fcff3a0df93895af759b15c90e3e).
