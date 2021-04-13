#!/usr/bin/env Rscript

# DM.R for DMRichR
# Author: Ben Laufer
# Contributors: Hyeyeon Hwang and Charles Mordaunt

# Initialize --------------------------------------------------------------

cat("\n[DMRichR] Initializing \t\t\t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")

options(scipen = 999)
options(readr.num_columns = 0)

if(length(grep("genomecenter.ucdavis.edu", .libPaths())) > 0){
  .libPaths("/share/lasallelab/programs/DMRichR/R_3.6")
  AnnotationHub::setAnnotationHubOption("CACHE", "/share/lasallelab/programs/DMRichR/R_3.6")
  ExperimentHub::setExperimentHubOption("CACHE", "/share/lasallelab/programs/DMRichR/R_3.6")
}else{
  sink("DMRichR_log.txt", type = "output", append = FALSE, split = TRUE)
}

if(suppressPackageStartupMessages(!require("DMRichR", quietly = TRUE))){
  Sys.setenv("R_REMOTES_NO_ERRORS_FROM_WARNINGS" = TRUE)
  BiocManager::install("ben-laufer/DMRichR")
  suppressPackageStartupMessages(library(DMRichR))
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
  optparse::make_option(c("-k", "--sexCheck"), type = "logical", default = FALSE,
              help = "Logical to confirm sex of each sample [default = %default]"),
  optparse::make_option(c("-e", "--ensembl"), type = "logical", default = FALSE,
              help = "Logical to select Ensembl transcript annotation database [default = %default]"),
  optparse::make_option(c("-f", "--GOfuncR"), type = "logical", default = TRUE,
              help = "Logical to run GOfuncR GO analysis [default = %default]")
  )
opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))

glue::glue("Assigning arguments to global variables...")
genome <- opt$genome
coverage <- opt$coverage
perGroup <- opt$perGroup
minCpGs <- opt$minCpGs
maxPerms <- opt$maxPerms
maxBlockPerms <- opt$maxBlockPerms
cutoff <- opt$cutoff
testCovariate <- opt$testCovariate
if(!is.null(opt$adjustCovariate)){
  adjustCovariate <- opt$adjustCovariate %>%
    strsplit(";") %>%
    unlist() %>%
    as.character()
}else if(is.null(opt$adjustCovariate)){
  adjustCovariate <- opt$adjustCovariate
}
matchCovariate <- opt$matchCovariate
cores <- opt$cores
cellComposition <-opt$cellComposition
sexCheck <-opt$sexCheck
EnsDb <- opt$ensembl
GOfuncR <- opt$GOfuncR

# Check for requirements
stopifnot(!is.null(genome))
stopifnot(!is.null(testCovariate))
stopifnot(coverage >= 1)

# Check for minimum number of cores
if(cores < 3){
  stop(glue::glue("You have selected {cores} cores, which is less than the 3 needed for DMRichR"))
}

# Check for more permutations than samples
nSamples <-  openxlsx::read.xlsx("sample_info.xlsx", colNames = TRUE) %>%
  nrow()

if(nSamples < maxPerms){
  print(glue::glue("Warning: You have requested {maxPerms} permutations for the DMR analysis, \\
                   which is more than the {nSamples} samples you have. \\
                   maxPerms will now be changed to {nSamples}."))
  maxPerms <- nSamples
}

if(nSamples < maxBlockPerms){
  print(glue::glue("Warning: You have requested {maxBlockPerms} permutations for the block analysis, \\
                   which is more than the {nSamples} samples you have. \\
                   maxBlockPerms will now be changed to {nSamples}."))
  maxBlockPerms <- nSamples
}

rm(nSamples)

# Print
glue::glue("genome = {genome}")
glue::glue("coverage = {coverage}")
glue::glue("perGroup = {perGroup}")
glue::glue("minCpGs = {minCpGs}")
glue::glue("maxPerms = {maxPerms}")
glue::glue("maxBlockPerms = {maxBlockPerms}")
glue::glue("cutoff = {cutoff}")
glue::glue("testCovariate = {testCovariate}")
glue::glue("adjustCovariate = {adjustCovariate}")
glue::glue("matchCovariate = {matchCovariate}")
glue::glue("cores = {cores}")
glue::glue("cellComposition = {cellComposition}")
glue::glue("sexCheck = {sexCheck}")
glue::glue("ensembl = {EnsDb}")
glue::glue("GOfuncR = {GOfuncR}")

# Setup annotation databases ----------------------------------------------

cat("\n[DMRichR] Selecting annotation databases \t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")

DMRichR::annotationDatabases(genome = genome,
                             EnsDb = EnsDb)

glue::glue("Saving Rdata...")
dir.create("RData")
settings_env <- ls(all = TRUE)
save(list = settings_env, file = "RData/settings.RData")
#load("RData/settings.RData")

# Load and process samples ------------------------------------------------

bs.filtered <- DMRichR::processBismark(files = list.files(path = getwd(),
                                                          pattern = "*.txt.gz"),
                                       meta = openxlsx::read.xlsx("sample_info.xlsx",
                                                                  colNames = TRUE) %>%
                                         dplyr::mutate_if(is.character, as.factor),
                                       testCovariate = testCovariate,
                                       adjustCovariate = adjustCovariate,
                                       matchCovariate = matchCovariate,
                                       coverage = coverage,
                                       cores = cores,
                                       perGroup = perGroup,
                                       sexCheck = sexCheck)

glue::glue("Saving Rdata...")
save(bs.filtered, file = "RData/bismark.RData")
#load("RData/bismark.RData")

glue::glue("Building annotations for plotting...")
if(is(TxDb, "TxDb")){
  annoTrack <- dmrseq::getAnnot(genome)
}else if(is(TxDb, "EnsDb")){
  annoTrack <- GenomicRanges::GRangesList(CpGs = DMRichR::getCpGs(genome),
                                          Exons = DMRichR::getExons(TxDb),
                                          compress = FALSE)
}

# Background --------------------------------------------------------------

cat("\n[DMRichR] Getting bsseq background regions \t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
dir.create("Extra")

DMRichR::getBackground(bs.filtered,
                       minNumRegion = minCpGs,
                       maxGap = 1000) %>% 
  write.table(file = "Extra/bsseq_background.csv",
              sep = ",",
              quote = FALSE,
              row.names = FALSE)

# Blocks ------------------------------------------------------------------

cat("\n[DMRichR] Testing for blocks with dmrseq \t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
start_time <- Sys.time()

tryCatch({
  blocks <- dmrseq::dmrseq(bs = bs.filtered,
                           cutoff = cutoff,
                           maxPerms = maxBlockPerms,
                           testCovariate = testCovariate,
                           adjustCovariate = adjustCovariate,
                           matchCovariate = matchCovariate,
                           block = TRUE,
                           minInSpan = 500,
                           bpSpan = 5e4,
                           maxGapSmooth = 1e6,
                           maxGap = 5e3,
                           minNumRegion = (minCpGs*2),
                           BPPARAM = BiocParallel::MulticoreParam(workers = cores)
  )
  
  glue::glue("Selecting significant blocks...")
  
  if(sum(blocks$qval < 0.05) == 0 & sum(blocks$pval < 0.05) != 0){
    sigBlocks <- blocks[blocks$pval < 0.05,]
  }else if(sum(blocks$qval < 0.05) >= 1){
    sigBlocks <- blocks[blocks$qval < 0.05,]
  }else if(sum(blocks$pval < 0.05) == 0 & length(blocks) != 0){
    glue::glue("No significant blocks detected in {length(blocks)} background blocks")
  }else if(length(blocks) == 0){
    glue::glue("No background blocks detected")
  }
  
  if(length(blocks) != 0){
    glue::glue("Exporting block and background information...")
    dir.create("Blocks")
    gr2csv(blocks, "Blocks/backgroundBlocks.csv")
    gr2bed(blocks, "Blocks/backgroundBlocks.bed")
    if(sum(blocks$pval < 0.05) > 0){
      glue::glue("{length(sigBlocks)} significant blocks of differential methylation in {length(blocks)} background blocks")
      gr2csv(sigBlocks, "Blocks/blocks.csv")
      gr2bed(sigBlocks, "Blocks/blocks.bed")
      
      glue::glue("Annotating and plotting blocks...")
      pdf("Blocks/Blocks.pdf", height = 7.50, width = 11.50)
      dmrseq::plotDMRs(bs.filtered,
                       regions = sigBlocks,
                       testCovariate = testCovariate,
                       annoTrack = annoTrack,
                       regionCol = "#FF00001A",
                       qval = FALSE,
                       stat = FALSE)
      dev.off()
    }
  }
  
  glue::glue("Blocks timing...")
  end_time <- Sys.time()
  end_time - start_time
  
  # Annotate blocks with gene symbols ---------------------------------------
  
  if(length(blocks) != 0){
    if(sum(blocks$pval < 0.05) > 0){
      glue::glue("Annotating blocks with gene symbols...")
      sigBlocks %>%
        DMRichR::annotateRegions(TxDb = TxDb,
                                 annoDb = annoDb) %T>%
        DMRichR::DMReport(regions = blocks,
                          bs.filtered = bs.filtered,
                          coverage = coverage,
                          name = "blockReport") %>% 
        openxlsx::write.xlsx(file = "Blocks/Blocks_annotated.xlsx")
    }
    
    glue::glue("Annotating background blocks with gene symbols...")
    blocks %>%
      DMRichR::annotateRegions(TxDb = TxDb,
                               annoDb = annoDb) %>% 
      openxlsx::write.xlsx(file = "Blocks/background_blocks_annotated.xlsx")
  }
  
  glue::glue("Saving RData...")
  save(blocks, file = "RData/Blocks.RData")
  #load("RData/Blocks.RData")
  
},
error = function(error_condition) {
  print(glue::glue("Block analysis has produced an error"))
})

# DMRs --------------------------------------------------------------------

cat("\n[DMRichR] Testing for DMRs with dmrseq \t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
start_time <- Sys.time()

regions <- dmrseq::dmrseq(bs = bs.filtered,
                          cutoff = cutoff,
                          minNumRegion = minCpGs,
                          maxPerms = maxPerms,
                          testCovariate = testCovariate,
                          adjustCovariate = adjustCovariate,
                          matchCovariate = matchCovariate,
                          BPPARAM = BiocParallel::MulticoreParam(workers = cores)
                          )

glue::glue("Selecting significant DMRs...", "\n")
if(sum(regions$qval < 0.05) < 100 & sum(regions$pval < 0.05) != 0){
  sigRegions <- regions[regions$pval < 0.05,]
}else if(sum(regions$qval < 0.05) >= 100){
  sigRegions <- regions[regions$qval < 0.05,]
}else if(sum(regions$pval < 0.05) == 0){
  stop(glue::glue("No significant DMRs detected in {length(regions)} background regions"))
}

glue::glue("Exporting DMR and background region information...")
dir.create("DMRs")
gr2bed(sigRegions, "DMRs/DMRs.bed")
gr2bed(regions, "DMRs/backgroundRegions.bed")

saveExternal(sigRegions = sigRegions,
             regions = regions)

if(sum(sigRegions$stat > 0) > 0 & sum(sigRegions$stat < 0) > 0){

  glue::glue("Summary: There are {tidySigRegions} DMRs \\
             ({tidyHyper}% hypermethylated, {tidyHypo}% hypomethylated) \\
             from {tidyRegions} background regions consisting of {tidyCpGs} CpGs \\
             assayed at {coverage}x coverage", 
             tidySigRegions = length(sigRegions),
             tidyHyper = round(sum(sigRegions$stat > 0) / length(sigRegions), digits = 2)*100,
             tidyHypo = round(sum(sigRegions$stat < 0) / length(sigRegions), digits = 2)*100,
             tidyRegions = length(regions),
             tidyCpGs = nrow(bs.filtered)
             )
}

glue::glue("DMR timing...")
end_time <- Sys.time()
end_time - start_time

glue::glue("Saving Rdata...")
save(regions, sigRegions, file = "RData/DMRs.RData")
#load("RData/DMRs.RData")

glue::glue("Annotating DMRs and plotting...")

pdf("DMRs/DMRs.pdf", height = 4, width = 8)
tryCatch({
  DMRichR::plotDMRs2(bs.filtered,
                     regions = sigRegions,
                     testCovariate = testCovariate,
                     extend = (end(sigRegions) - start(sigRegions) + 1)*2,
                     addRegions = sigRegions,
                     annoTrack = annoTrack,
                     regionCol = "#FF00001A",
                     lwd = 2,
                     qval = FALSE,
                     stat = FALSE,
                     horizLegend = FALSE)
  },
  error = function(error_condition) {
    print(glue::glue("ERROR: One (or more) of your DMRs can't be plotted, \\
                      try again later by manually loading R Data and subsetting sigRegions"))
  })
dev.off()

# Annotate DMRs with gene symbols -----------------------------------------

cat("\n[DMRichR] Annotating DMRs with gene symbols \t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")

sigRegions %>%
  DMRichR::annotateRegions(TxDb = TxDb,
                           annoDb = annoDb) %T>%
  DMRichR::DMReport(regions = regions,
                    bs.filtered = bs.filtered,
                    coverage = coverage,
                    name = "DMReport") %>% 
  openxlsx::write.xlsx(file = "DMRs/DMRs_annotated.xlsx")

glue::glue("Annotating background regions with gene symbols...")
regions %>%
  DMRichR::annotateRegions(TxDb = TxDb,
                           annoDb = annoDb) %>% 
  openxlsx::write.xlsx(file = "DMRs/background_annotated.xlsx")

# Individual smoothed values ----------------------------------------------

cat("\n[DMRichR] Smoothing individual methylation values \t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
start_time <- Sys.time()

bs.filtered.bsseq <- bsseq::BSmooth(bs.filtered,
                                    BPPARAM = BiocParallel::MulticoreParam(workers = cores,
                                                                           progressbar = TRUE)
                                    )

# Drop chrY in Rat only due to poor quality (some CpGs in females map to Y)
if(genome == "rn6"){
  bs.filtered.bsseq <- GenomeInfoDb::dropSeqlevels(bs.filtered.bsseq,
                                                   "chrY",
                                                   pruning.mode = "coarse")
  GenomeInfoDb::seqlevels(bs.filtered.bsseq)
}

bs.filtered.bsseq

glue::glue("Extracting individual smoothed methylation values of DMRs...")
bs.filtered.bsseq %>%
  DMRichR::getSmooth(regions) %>%
  DMRichR::smooth2txt("DMRs/DMR_individual_smoothed_methylation.txt")

glue::glue("Extracting individual smoothed methylation values of background regions for WGCNA...")
bs.filtered.bsseq %>%
  DMRichR::getSmooth(regions) %>%
  DMRichR::smooth2txt("DMRs/background_region_individual_smoothed_methylation.txt")

glue::glue("Individual smoothing timing...")
end_time <- Sys.time()
end_time - start_time

glue::glue("Saving Rdata...")
save(bs.filtered.bsseq,
     file = "RData/bsseq.RData")
#load("RData/bsseq.RData")

# ChromHMM and Roadmap Epigenomics ----------------------------------------

if(length(grep("genomecenter.ucdavis.edu", .libPaths())) > 0 & genome == "hg38"){
  
  dir.create("Extra/LOLA")
  setwd("Extra/LOLA")
  
  dmrList <- sigRegions %>% 
    DMRichR::dmrList()
  
  LOLA <- function(x){
    
    dir.create(names(dmrList)[x])
    setwd(names(dmrList)[x])
    
    dmrList[x] %>%
      DMRichR::chromHMM(regions = regions,
                        cores = floor(cores/3)) %>% 
      DMRichR::chromHMM_heatmap()
    
    dmrList[x] %>%
      DMRichR::roadmap(regions = regions,
                       cores = floor(cores/3)) %>% 
      DMRichR::roadmap_heatmap()
    
    if(file.exists("Rplots.pdf")){file.remove("Rplots.pdf")}
  }
  
  parallel::mclapply(seq_along(dmrList),
                     LOLA,
                     mc.cores = 3,
                     mc.silent = TRUE)
  
  setwd("../..")
}

# Smoothed global, chromosomal, and CGi methylation statistics ------------

dir.create("Global")

bs.filtered.bsseq %>%
  DMRichR::globalStats(genome = genome,
                       testCovariate = testCovariate,
                       adjustCovariate = adjustCovariate,
                       matchCovariate = matchCovariate) %>%
  openxlsx::write.xlsx("Global/smoothed_globalStats.xlsx") 

# Global plots ------------------------------------------------------------

windows <- bs.filtered.bsseq %>%
  DMRichR::windows(goi = goi)

CpGs <- bs.filtered.bsseq %>%
  DMRichR::CpGs()

plots <- c("windows", "CpGs")

if(genome %in% c("hg38", "hg19", "mm10", "mm9", "rheMac10", "rheMac8", "rn6", "danRer11", "galGal6",
                 "bosTau9", "panTro6", "dm6", "susScr11", "canFam3")){
  
  CGi <- bs.filtered.bsseq %>% 
    DMRichR::CGi(genome = genome)
  
  plots <- c("windows", "CpGs", "CGi")
}


purrr::walk(plots,
            function(plotMatrix,
                     group =  bs.filtered.bsseq %>%
                       pData() %>%
                       dplyr::as_tibble() %>%
                       dplyr::pull(!!testCovariate) %>%
                       forcats::fct_rev()){
              
              title <- dplyr::case_when(plotMatrix == "windows" ~ "20Kb Windows",
                                        plotMatrix == "CpGs" ~ "Single CpG",
                                        plotMatrix == "CGi" ~ "CpG Island")
              
              plotMatrix %>%
                get() %>% 
                DMRichR::PCA(group = group) %>%
                ggplot2::ggsave(glue::glue("Global/{title} PCA.pdf"),
                                plot = .,
                                device = NULL,
                                width = 11,
                                height = 8.5)
              
              plotMatrix %>%
                get() %>% 
                DMRichR::densityPlot(group = group) %>% 
                ggplot2::ggsave(glue::glue("Global/{title} Density Plot.pdf"),
                                plot = .,
                                device = NULL,
                                width = 11,
                                height = 4)
              
              Glimma::glMDSPlot(plotMatrix %>%
                                  get(),
                                groups = cbind(bsseq::sampleNames(bs.filtered.bsseq),
                                               pData(bs.filtered.bsseq)) %>%
                                  dplyr::as_tibble() %>% 
                                  dplyr::select(-col) %>%
                                  dplyr::rename(Name = bsseq..sampleNames.bs.filtered.bsseq.),
                                path = getwd(),
                                folder = "interactiveMDS",
                                html = glue::glue("{title} MDS plot"),
                                launch = FALSE)
            })

# Heatmap -----------------------------------------------------------------

sigRegions %>%
  DMRichR::smoothPheatmap(bs.filtered.bsseq = bs.filtered.bsseq,
                          testCovariate = testCovariate)

# CpG and genic enrichment testing ----------------------------------------

cat("\n[DMRichR] Performing DMRichments \t\t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")

DMRich <- function(x){
  
  if(genome %in% c("hg38", "hg19", "mm10", "mm9", "rheMac10", "rheMac8", "rn6", "danRer11", "galGal6", "bosTau9", "panTro6", "dm6", "susScr11", "canFam3")){
    print(glue::glue("Running CpG annotation enrichments for {names(dmrList)[x]}"))
    dmrList[x] %>% 
      DMRichR::DMRichCpG(regions = regions,
                         genome = genome) %T>%
      openxlsx::write.xlsx(file = glue::glue("DMRichments/{names(dmrList)[x]}_CpG_enrichments.xlsx")) %>% 
      DMRichR::DMRichPlot(type = "CpG") %>% 
      ggplot2::ggsave(glue::glue("DMRichments/{names(dmrList)[x]}_CpG_enrichments.pdf"),
                      plot = ., 
                      width = 4,
                      height = 3)
  }
  
  print(glue::glue("Running gene region annotation enrichments for {names(dmrList)[x]}"))
  dmrList[x] %>% 
    DMRichR::DMRichGenic(regions = regions,
                         TxDb = TxDb,
                         annoDb = annoDb) %T>%
    openxlsx::write.xlsx(file = glue::glue("DMRichments/{names(dmrList)[x]}_genic_enrichments.xlsx")) %>% 
    DMRichR::DMRichPlot(type = "genic") %>% 
    ggplot2::ggsave(glue::glue("DMRichments/{names(dmrList)[x]}_genic_enrichments.pdf"),
                    plot = ., 
                    width = 4,
                    height = 4)
}

dmrList <- sigRegions %>% 
  DMRichR::dmrList()

dir.create("DMRichments")

purrr::walk(seq_along(dmrList),
            DMRich)

purrr::walk(dplyr::case_when(genome %in% c("hg38", "hg19", "mm10", "mm9", "rn6") ~ c("CpG", "genic"),
                             TRUE ~ "genic") %>%
              unique(),
            function(type){
              
              print(glue::glue("Creating DMRichMultiPlots for {type} annotations"))
              
              DMRichR::DMparseR(direction =  c("All DMRs",
                                               "Hypermethylated DMRs",
                                               "Hypomethylated DMRs"),
                                type = type) %>%
                DMRichR::DMRichPlot(type = type,
                                    multi = TRUE) %>% 
                ggplot2::ggsave(glue::glue("DMRichments/{type}_multi_plot.pdf"),
                                plot = .,
                                device = NULL,
                                height = dplyr::case_when(type == "genic" ~ 5,
                                                          type == "CpG" ~ 3.5),
                                width = 7)
            })

# Overlap with human imprinted genes --------------------------------------

cat("\n[DMRichR] Testing for imprinted gene enrichment \t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")

dmrList <- sigRegions %>% 
  DMRichR::dmrList()

sink("DMRs/human_imprinted_gene_overlaps.txt")

purrr::walk(seq_along(dmrList),
            function(x){
              print(glue::glue("Analyzing {names(dmrList)[x]}"))
              
              dmrList[x] %>%
                DMRichR::imprintOverlap(regions = regions,
                                        TxDb = TxDb,
                                        annoDb = annoDb)
            })

sink()

# Manhattan and Q-Q plots -------------------------------------------------

regions %>%
  DMRichR::annotateRegions(TxDb = TxDb,
                           annoDb = annoDb) %>% 
  DMRichR::manQQ()

# Gene Ontology analyses --------------------------------------------------

cat("\n[DMRichR] Performing gene ontology analyses \t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")

dir.create("Ontologies")

if(genome %in% c("hg38", "hg19", "mm10", "mm9")){
  
  print(glue::glue("Running GREAT"))
  GREATjob <- sigRegions %>%
    dplyr::as_tibble() %>%
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE) %>%
    rGREAT::submitGreatJob(bg = regions,
                           species = genome,
                           rule = "oneClosest",
                           request_interval = 1,
                           version = "4.0.4")
  
  print(glue::glue("Saving and plotting GREAT results"))
  GREATjob %>%
    rGREAT::getEnrichmentTables(category = "GO") %T>% #%>%
    #purrr::map(~ dplyr::filter(., Hyper_Adjp_BH < 0.05)) %T>%
    openxlsx::write.xlsx(file = glue::glue("Ontologies/GREAT_results.xlsx")) %>%
    DMRichR::slimGO(tool = "rGREAT",
                    annoDb = annoDb,
                    plots = FALSE) %T>%
    openxlsx::write.xlsx(file = glue::glue("Ontologies/GREAT_slimmed_results.xlsx")) %>%
    DMRichR::GOplot() %>%
    ggplot2::ggsave(glue::glue("Ontologies/GREAT_plot.pdf"),
                    plot = .,
                    device = NULL,
                    height = 12,
                    width = 16)
  
  # pdf(glue::glue("Ontologies/GREAT_gene_associations_graph.pdf"),
  #     height = 8.5,
  #     width = 11)
  # par(mfrow = c(1, 3))
  # res <- rGREAT::plotRegionGeneAssociationGraphs(GREATjob)
  # dev.off()
  # write.csv(as.data.frame(res),
  #           file = glue::glue("Ontologies/GREATannotations.csv"),
  #           row.names = FALSE)
}

if(GOfuncR == TRUE){
  print(glue::glue("Running GOfuncR"))
  sigRegions %>% 
    DMRichR::GOfuncR(regions = regions,
                     n_randsets = 1000,
                     upstream = 5000,
                     downstream = 1000,
                     annoDb = annoDb,
                     TxDb = TxDb) %T>%
    openxlsx::write.xlsx(glue::glue("Ontologies/GOfuncR.xlsx")) %>% 
    DMRichR::slimGO(tool = "GOfuncR",
                    annoDb = annoDb,
                    plots = FALSE) %T>%
    openxlsx::write.xlsx(file = glue::glue("Ontologies/GOfuncR_slimmed_results.xlsx")) %>% 
    DMRichR::GOplot() %>% 
    ggplot2::ggsave(glue::glue("Ontologies/GOfuncR_plot.pdf"),
                    plot = .,
                    device = NULL,
                    height = 12,
                    width = 16)
}

if(genome != "TAIR10" & genome != "TAIR9"){
  print(glue::glue("Running enrichR"))
  
  enrichR:::.onAttach() # Needed or else "EnrichR website not responding"
  #dbs <- enrichR::listEnrichrDbs()
  dbs <- c("GO_Biological_Process_2018",
           "GO_Cellular_Component_2018",
           "GO_Molecular_Function_2018",
           "KEGG_2019_Human",
           "Panther_2016",
           "Reactome_2016",
           "RNA-Seq_Disease_Gene_and_Drug_Signatures_from_GEO")
  
  if(genome %in% c("mm10", "mm9", "rn6")){
    dbs %>%
      gsub(pattern = "Human", replacement = "Mouse")
  }else if(genome %in% c("danRer11", "dm6")){
    if(genome == "danRer11"){
      enrichR::setEnrichrSite("FishEnrichr")
    }else if(genome == "dm6"){
      enrichR::setEnrichrSite("FlyEnrichr")}
    dbs <- c("GO_Biological_Process_2018",
             "GO_Cellular_Component_2018",
             "GO_Molecular_Function_2018",
             "KEGG_2019")
  }
  
  sigRegions %>%
    DMRichR::annotateRegions(TxDb = TxDb,
                             annoDb = annoDb) %>%  
    dplyr::select(geneSymbol) %>%
    purrr::flatten() %>%
    enrichR::enrichr(dbs) %T>% #%>% 
    #purrr::map(~ dplyr::filter(., Adjusted.P.value < 0.05)) %T>%
    openxlsx::write.xlsx(file = glue::glue("Ontologies/enrichr.xlsx")) %>%
    DMRichR::slimGO(tool = "enrichR",
                    annoDb = annoDb,
                    plots = FALSE) %T>%
    openxlsx::write.xlsx(file = glue::glue("Ontologies/enrichr_slimmed_results.xlsx")) %>% 
    DMRichR::GOplot() %>% 
    ggplot2::ggsave(glue::glue("Ontologies/enrichr_plot.pdf"),
                    plot = .,
                    device = NULL,
                    height = 12,
                    width = 16)
}

# Machine learning --------------------------------------------------------

methylLearnOutput <- DMRichR::methylLearn(bs.filtered.bsseq = bs.filtered.bsseq,
                                          sigRegions = sigRegions,
                                          testCovariate = testCovariate,
                                          TxDb = TxDb,
                                          annoDb = annoDb,
                                          topPercent = 1,
                                          output = "all",
                                          saveHtmlReport = TRUE)

if(!dir.exists("./Machine_learning")) {
  dir.create("./Machine_learning")
} 

if(length(methylLearnOutput) == 1) {
  openxlsx::write.xlsx(list(Annotations_Common_DMRs = methylLearnOutput), 
                       file = "./Machine_learning/Machine_learning_output_one.xlsx") 
} else {
  openxlsx::write.xlsx(list(Annotations_Common_DMRs = methylLearnOutput$`Annotated common DMRs`,
                            RF_Ranking_All_DMRs = methylLearnOutput$`RF ranking`,
                            SVM_Ranking_All_DMRs = methylLearnOutput$`SVM ranking`),
                       file = "./Machine_learning/Machine_learning_output_all.xlsx") 
}

glue::glue("Saving RData...")
save(methylLearnOutput, file = "RData/machineLearning.RData")
#load("RData/machineLearing.RData")

# End ---------------------------------------------------------------------

cat("\n[DMRichR] Summary \t\t\t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")

glue::glue("Summary: There are {tidySigRegions} \\
             ({tidyHyper}% hypermethylated, {tidyHypo}% hypomethylated) \\
             from {tidyRegions} background regions consisting of {tidyCpGs} CpGs \\
             assayed at {coverage}x coverage", 
           tidySigRegions = length(sigRegions),
           tidyHyper = round(sum(sigRegions$stat > 0) / length(sigRegions),
                             digits = 2)*100,
           tidyHypo = round(sum(sigRegions$stat < 0) / length(sigRegions),
                            digits = 2)*100,
           tidyRegions = length(regions),
           tidyCpGs = nrow(bs.filtered)
           )

if(sum(blocks$pval < 0.05) > 0 & length(blocks) != 0){
glue::glue("{length(sigBlocks)} significant blocks of differential methylation \\
           in {length(blocks)} background blocks")
}

if(length(grep("genomecenter.ucdavis.edu", .libPaths())) == 0){sessionInfo()}
if(file.exists("Rplots.pdf")){file.remove("Rplots.pdf")}
rm(list = ls())
glue::glue("Done...")
