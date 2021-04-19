#!/usr/bin/env Rscript

# DM.R for DMRichR
# Author: Ben Laufer
# Contributors: Hyeyeon Hwang and Charles Mordaunt

# Initialize --------------------------------------------------------------

cat("\n[DMRichR] Initializing \t\t\t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")

if(!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if(!requireNamespace("remotes", quietly = TRUE))
  install.packages("remotes")
if(suppressPackageStartupMessages(!requireNamespace("DMRichR", quietly = TRUE))){
  Sys.setenv("R_REMOTES_NO_ERRORS_FROM_WARNINGS" = TRUE)
  BiocManager::install("ben-laufer/DMRichR")
}
suppressPackageStartupMessages(library(DMRichR))

if(length(grep("genomecenter.ucdavis.edu", .libPaths())) > 0){
  .libPaths("/share/lasallelab/programs/DMRichR/R_3.6")
  AnnotationHub::setAnnotationHubOption("CACHE", "/share/lasallelab/programs/DMRichR/R_3.6")
  ExperimentHub::setExperimentHubOption("CACHE", "/share/lasallelab/programs/DMRichR/R_3.6")
}

# Global variables --------------------------------------------------------

cat("\n[DMRichR] Processing arguments from script \t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")

option_list <- list(
  optparse::make_option(c("-g", "--genome"), type = "character", default = NULL,
                        help = "Choose a genome (hg38, hg19, mm10, mm9, rheMac10, rheMac8, rn6, danRer11, galGal6, bosTau9, panTro6, dm6, susScr11, canFam3, TAIR10, or TAIR9) [required]"),
  optparse::make_option(c("-x", "--coverage"), type = "integer", default = 1,
                        help = "Choose a CpG coverage cutoff [default = %default]"),
  optparse::make_option(c("-s", "--perGroup"), type = "double", default = 1,
                        help = "Choose the percent [values from 0 to 1] of samples in all combinations of covariates meeting CpG coverage cutoff [default = %default]"),
  optparse::make_option(c("-n", "--minCpGs"), type = "integer", default = 5,
                        help = "Choose the minimum number of CpGs for a DMR [default = %default]"),
  optparse::make_option(c("-p", "--maxPerms"), type = "integer", default = 10,
                        help = "Choose the number of permutations for the DMR analysis [default = %default]"),
  optparse::make_option(c("-b", "--maxBlockPerms"), type = "integer", default = 10,
                        help = "Choose the number of permutations for the block analysis [default = %default]"),
  optparse::make_option(c("-o", "--cutoff"), type = "double", default = 0.05,
                        help = "Choose the cutoff value [from 0 to 1] for the single CpG coefficient utilized to discover testable background regions [default = %default]"),
  optparse::make_option(c("-t", "--testCovariate"), type = "character", default = NULL,
                        help = "Choose a test covariate [required]"),
  optparse::make_option(c("-a", "--adjustCovariate"), type = "character", default = NULL,
                        help = "Choose covariates to directly adjust [default = NULL]"),
  optparse::make_option(c("-m", "--matchCovariate"), type = "character", default = NULL,
                        help = "Choose covariate to balance permutations [default = NULL]"),
  optparse::make_option(c("-c", "--cores"), type = "integer", default = 20,
                        help = "Choose number of cores [default = %default]"),
  optparse::make_option(c("-e", "--cellComposition"), type = "logical", default = FALSE,
                        help = "Logical to estimate blood cell composition [default = %default]"),
  optparse::make_option(c("-k", "--sexCheck"), type = "logical", default = FALSE,
                        help = "Logical to confirm sex of each sample [default = %default]"),
  optparse::make_option(c("-d", "--EnsDb"), type = "logical", default = FALSE,
                        help = "Logical to select Ensembl transcript annotation database [default = %default]"),
  optparse::make_option(c("-f", "--GOfuncR"), type = "logical", default = TRUE,
                        help = "Logical to run GOfuncR GO analysis [default = %default]")
)
opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))

# DM.R --------------------------------------------------------------------

DMRichR::DM.R(genome = opt$genome,
              coverage = opt$coverage,
              perGroup = opt$perGroup,
              minCpGs =  opt$minCpGs,
              maxPerms =  opt$maxPerms,
              maxBlockPerms = opt$maxBlockPerms,
              cutoff = opt$cutoff,
              testCovariate = opt$testCovariate,
              adjustCovariate = if(!is.null(opt$adjustCovariate)){
                adjustCovariate <- opt$adjustCovariate %>%
                  strsplit(";") %>%
                  unlist() %>%
                  as.character()
                }else if(is.null(opt$adjustCovariate)){
                  adjustCovariate <- opt$adjustCovariate
                  },
              matchCovariate = opt$matchCovariate,
              cores = opt$cores,
              GOfuncR = opt$GOfuncR,
              sexCheck = opt$sexCheck,
              EnsDb = opt$EnsDb,
              cellComposition = opt$cellComposition)
    