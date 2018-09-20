# DM.R
A workflow for the statistical analysis and visualization of differentially methylated regions (DMRs) from a CpG count matrix

Use adjustCovariate for covariates that are continuous or contain two or more groups. More than one covariate can be adjusted for.
Use matchCovariate for balancing permutations, which is ideal for two group covariates such as sex. Only one covariate can be balanced.

Keegan re beta coefficent: You’re exactly right that it represents the average effect size over the region, but if you’d like to take it a step further and connect it to the difference seen in the plot, you can divide the beta coefficient by pi (yep, 3.14159…) to put it on the scale of a proportion difference. This is because the beta coefficient is on the scale of the arcsine transformed differences. So beta/pi will be similar to (and correlated with) the simple mean proportion difference across the region, but the beta/pi quantity from the model is adjusted for things like coverage and correlated errors. 

Keegan re individual values: The per sample smoothing lines in the plots (1) are very different than the smoothed methylation differences dmrseq computes,and (2) are *purely* for visualization purposes. They simply smooth the methylation values with loess, and do not use the model in any way. If you really need smoothed sample-specific methylation values, I’d suggest obtaining them with the bsseq package.
