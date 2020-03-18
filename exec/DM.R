#!/usr/bin/env Rscript

# DMRichR
# Author: Ben Laufer
# Contributors: Charles Mordaunt and Hyeyeon Hwang

# R settings --------------------------------------------------------------

rm(list=ls())
options(scipen=999)

if(length(grep("genomecenter.ucdavis.edu", .libPaths())) > 0){
  .libPaths("/share/lasallelab/programs/DMRichR/R_3.6")
  AnnotationHub::setAnnotationHubOption("CACHE", "/share/lasallelab/programs/DMRichR/R_3.6")
}else{
  sink("DMRichR_log.txt", type = "output", append = FALSE, split = TRUE)
}

# Install and update ------------------------------------------------------

cat("\n[DMRichR] Installing and updating packages \t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")

if (!requireNamespace(c("BiocManager", "remotes"), quietly = TRUE))
  install.packages(c("BiocManager", "remotes"))
suppressWarnings(BiocManager::valid(fix = TRUE, update = TRUE, ask = FALSE))
BiocManager::install("ben-laufer/DMRichR", update = TRUE, ask = FALSE)
suppressPackageStartupMessages(library(DMRichR, attach.required = T))

# Global variables --------------------------------------------------------

cat("\n[DMRichR] Processing arguments from script \t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")

option_list <- list(
  optparse::make_option(c("-g", "--genome"), type = "character", default = NULL,
              help = "Choose a genome (hg38, hg19, mm10, mm9, rheMac10, rheMac8, rn6) [required]"),
  optparse::make_option(c("-x", "--coverage"), type = "integer", default = 1,
              help = "Choose a CpG coverage cutoff [default = %default]"),
  optparse::make_option(c("-s", "--perGroup"), type = "double", default = 1,
              help = "Choose the percent [values from 0 to 1] of samples in all combinations of covariates meeting CpG coverage cutoff [default = %default]"),
  optparse::make_option(c("-n", "--minCpGs"), type = "integer", default = 5,
              help = "Choose the minimum number of CpGs for a DMR [default = %default]"),
  optparse::make_option(c("-p", "--maxPerms"), type = "integer", default = 10,
              help = "Choose the number of permutations for DMR and block analyses [default = %default]"),
  optparse::make_option(c("-o", "--cutoff"), type = "double", default = 0.05,
              help = "Choose the cutoff value [from 0 to 1] for the single CpG coefficient utilized to discover testable background regions [default = %default]"),
  optparse::make_option(c("-t", "--testCovariate"), type = "character", default = NULL,
              help = "Choose a test covariate [required]"),
  optparse::make_option(c("-a", "--adjustCovariate"), type = "character", default = NULL,
              help = "Choose covariates to directly adjust [default = NULL]"),
  optparse::make_option(c("-m", "--matchCovariate"), type = "character", default = NULL,
              help = "Choose covariate to balance permutations [default = NULL]"),
  optparse::make_option(c("-c", "--cores"), type = "integer", default = 20,
              help = "Choose number of cores [default = %default]")
)
opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))

glue::glue("Assigning arguments to global variables...")
# Check for requirements
stopifnot(!is.null(opt$genome))
stopifnot(!is.null(opt$testCovariate))
stopifnot(opt$coverage >= 1)
# Assign
genome <- as.character(opt$genome)
coverage <- as.numeric(opt$coverage)
perGroup <- as.numeric(opt$perGroup)
minCpGs <- as.numeric(opt$minCpGs)
maxPerms <- as.numeric(opt$maxPerms)
cutoff <- as.numeric(opt$cutoff)
testCovariate <- as.character(opt$testCovariate)
if(!is.null(opt$adjustCovariate)){
  adjustCovariate <- opt$adjustCovariate %>%
    strsplit(";") %>%
    unlist() %>%
    as.character()
}else if(is.null(opt$adjustCovariate)){
  adjustCovariate <- opt$adjustCovariate
}
if(!is.null(opt$matchCovariate)){
  matchCovariate <- as.character(opt$matchCovariate)
}else if(is.null(opt$matchCovariate)){
  matchCovariate <- opt$matchCovariate
}
cores <- as.numeric(opt$cores)

# Print
glue::glue("genome = {genome}")
glue::glue("coverage = {coverage}")
glue::glue("perGroup = {perGroup}")
glue::glue("minCpGs = {minCpGs}")
glue::glue("maxPerms = {maxPerms}")
glue::glue("cutoff = {cutoff}")
glue::glue("testCovariate = {testCovariate}")
glue::glue("adjustCovariate = {adjustCovariate}")
glue::glue("matchCovariate = {matchCovariate}")
glue::glue("cores = {cores}")

# Setup annotation databases ----------------------------------------------

cat("\n[DMRichR] Selecting annotation databases \t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")

packages <- dplyr::case_when(genome == "hg38" ~ c("BSgenome.Hsapiens.UCSC.hg38", "TxDb.Hsapiens.UCSC.hg38.knownGene", "org.Hs.eg.db"),
                             genome == "hg19" ~ c("BSgenome.Hsapiens.UCSC.hg19", "TxDb.Hsapiens.UCSC.hg19.knownGene", "org.Hs.eg.db"),
                             genome == "mm10" ~ c("BSgenome.Mmusculus.UCSC.mm10", "TxDb.Mmusculus.UCSC.mm10.knownGene", "org.Mm.eg.db"),
                             genome == "mm9" ~ c("BSgenome.Mmusculus.UCSC.mm9", "TxDb.Mmusculus.UCSC.mm9.knownGene", "org.Mm.eg.db"),
                             genome == "rheMac10" ~ c("BSgenome.Mmulatta.UCSC.rheMac10", "TxDb.Mmulatta.UCSC.rheMac10.refGene", "org.Mmu.eg.db"),
                             genome == "rheMac8" ~ c("BSgenome.Mmulatta.UCSC.rheMac8", "TxDb.Mmulatta.UCSC.rheMac8.refGene", "org.Mmu.eg.db"),
                             genome == "rn6" ~ c("BSgenome.Rnorvegicus.UCSC.rn6", "TxDb.Rnorvegicus.UCSC.rn6.refGene", "org.Rn.eg.db")
                             )

new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new.packages)>0){
  glue::glue("Installing {packages}")
  suppressMessages(BiocManager::install(new.packages, ask = FALSE, quiet = TRUE))
  cat("Done", "\n")
}
glue::glue("Loading {packages}")
stopifnot(suppressMessages(sapply(packages, require, character.only = TRUE)))

if(genome == "hg38"){
  goi <- BSgenome.Hsapiens.UCSC.hg38
  TxDb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  annoDb <- "org.Hs.eg.db"
}else if(genome == "hg19"){
  goi <- BSgenome.Hsapiens.UCSC.hg19 
  TxDb <- TxDb.Hsapiens.UCSC.hg19.knownGene
  annoDb <- "org.Hs.eg.db"
}else if(genome == "mm10"){
  goi <- BSgenome.Mmusculus.UCSC.mm10
  TxDb <- TxDb.Mmusculus.UCSC.mm10.knownGene
  annoDb <- "org.Mm.eg.db"
}else if(genome == "mm9"){
  goi <- BSgenome.Mmusculus.UCSC.mm9
  TxDb <- TxDb.Mmusculus.UCSC.mm9.knownGene
  annoDb <- "org.Mm.eg.db"
}else if(genome == "rheMac10"){
  goi <- BSgenome.Mmulatta.UCSC.rheMac10
  TxDb <- TxDb.Mmulatta.UCSC.rheMac10.refGene
  annoDb <- "org.Mmu.eg.db"
}else if(genome == "rheMac8"){
  goi <- BSgenome.Mmulatta.UCSC.rheMac8
  TxDb <- TxDb.Mmulatta.UCSC.rheMac8.refGene
  annoDb <- "org.Mmu.eg.db"
}else if(genome == "rn6"){
  goi <- BSgenome.Rnorvegicus.UCSC.rn6
  TxDb <- TxDb.Rnorvegicus.UCSC.rn6.refGene
  annoDb <- "org.Rn.eg.db"
}else{
  stop(glue("{genome} is not supported, please choose either hg38, hg19, mm10, mm9, rheMac10, rheMac8, or rn6 [Case Sensitive]"))
}

# Load and process samples ------------------------------------------------

bs.filtered <- processBismark(files = list.files(path = getwd(), pattern = "*.txt.gz"),
                              meta = openxlsx::read.xlsx("sample_info.xlsx", colNames = TRUE) %>%
                                dplyr::mutate_if(is.character, as.factor),
                              testCovar = testCovariate,
                              adjustCovar = adjustCovariate,
                              matchCovar = matchCovariate,
                              Cov = coverage,
                              mc.cores = cores,
                              per.Group = perGroup)

glue::glue("Saving Rdata...")
dir.create("RData")
bismark_env <- ls(all = TRUE)
save(list = bismark_env, file = "RData/bismark.RData")
#load("RData/bismark.RData")

# Distribution plots ------------------------------------------------------

#cat("\n[DMRichR] Plotting Empirical Distribution of CpGs \t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
#pdf("Filtered_CpG_Methylation_Distributions.pdf", height = 7.50, width = 11.50)
#plotEmpiricalDistribution(bs.filtered, testCovariate = testCovariate)
#plotEmpiricalDistribution(bs.filtered, testCovariate = testCovariate, type = "Cov", bySample = TRUE)
#dev.off()

# Background --------------------------------------------------------------

cat("\n[DMRichR] Getting bsseq background regions \t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
dir.create("Extra")

getBackground(bs.filtered,
              minNumRegion = minCpGs,
              maxGap = 1000) %>% 
  write.table(file = "Extra/bsseq_background.csv",
              sep = ",",
              quote = FALSE,
              row.names = FALSE)

# Blocks ------------------------------------------------------------------

cat("\n[DMRichR] Testing for blocks of differential methylation", format(Sys.time(), "%d-%m-%Y %X"), "\n")
start_time <- Sys.time()

blocks <- dmrseq::dmrseq(bs = bs.filtered,
                         cutoff = cutoff,
                         maxPerms = (maxPerms*3),
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
                     annoTrack = getAnnot(genome),
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
      annotateRegions(TxDb = TxDb,
                      annoDb = annoDb) %T>%
      DMReport(regions = blocks,
               bsseq = bs.filtered.bsseq,
               coverage = coverage,
               name = "blockReport") %>% 
      openxlsx::write.xlsx(file = "Blocks/Blocks_annotated.xlsx")
  }
  
  glue::glue("Annotating background blocks with gene symbols...")
  blocks %>%
    annotateRegions(TxDb = TxDb,
                    annoDb = annoDb) %>% 
    openxlsx::write.xlsx(file = "Blocks/background_blocks_annotated.xlsx")
}

glue::glue("Saving RData...")
blocks_env <- ls(all = TRUE)[!(ls(all = TRUE) %in% bismark_env)]
save(list = blocks_env, file = "RData/Blocks.RData")
#load("RData/Blocks.RData")

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
  glue::glue("{length(sigRegions)} Significant DMRs \\
             ({round(sum(sigRegions$stat > 0) / length(sigRegions), digits = 2)*100}% hypermethylated, \\
             {round(sum(sigRegions$stat < 0) / length(sigRegions), digits = 2)*100}% hypomethylated) \\
             in {length(regions)} background regions \\
             from {nrow(bs.filtered)} CpGs assayed at {coverage}x coverage")
}

glue::glue("Saving Rdata...")
DMRs_env <- ls(all = TRUE)[!(ls(all = TRUE) %in% bismark_env) &
                             !(ls(all = TRUE) %in% blocks_env)]
save(list = DMRs_env, file = "RData/DMRs.RData")
#load("RData/DMRs.RData")

glue::glue("DMR timing...")
end_time <- Sys.time()
end_time - start_time

# Individual smoothed values ----------------------------------------------

cat("\n[DMRichR] Smoothing individual methylation values \t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
start_time <- Sys.time()

bs.filtered.bsseq <- BSmooth(bs.filtered,
                             BPPARAM = MulticoreParam(workers = cores,
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
  getSmooth(regions) %>%
  smooth2txt("DMRs/DMR_individual_smoothed_methylation.txt")

glue::glue("Extracting individual smoothed methylation values of background regions for WGCNA...")
bs.filtered.bsseq %>%
  getSmooth(regions) %>%
  smooth2txt("DMRs/background_region_individual_smoothed_methylation.txt")

glue::glue("Saving Rdata...")
bsseq_env <- ls(all = TRUE)[!(ls(all = TRUE) %in% bismark_env) &
                              !(ls(all = TRUE) %in% blocks_env) &
                              !(ls(all = TRUE) %in% DMRs_env)]
save(list = bsseq_env, file = "RData/bsseq.RData")
#load("RData/bsseq.RData")

glue::glue("Individual smoothing timing...")
end_time <- Sys.time()
end_time - start_time

# Plot smoothed DMR methylation -------------------------------------------

# glue::glue("Annotating and plotting...")
# pdf("DMRs.pdf", height = 7.50, width = 11.50)
# annoTrack <- getAnnot(genome)
# plotDMRs(bs.filtered,
#          regions = sigRegions,
#          testCovariate = testCovariate,
#          annoTrack = annoTrack,
#          qval = F)
# dev.off()

glue::glue("Annotating DMRs and plotting smoothed values...")
pData <- pData(bs.filtered.bsseq)
if(length(levels(pData[,testCovariate])) == 2){
  pData$col <- NULL
  pData$col[pData[,testCovariate] == levels(pData[,testCovariate])[1]] <- "#3366CC"
  pData$col[pData[,testCovariate] == levels(pData[,testCovariate])[2]] <- "#FF3366"
  pData(bs.filtered.bsseq) <- pData
}
 
pdf("DMRs/DMRs.pdf", height = 4, width = 8)
plotDMRs2(bs.filtered.bsseq,
          regions = sigRegions,
          testCovariate = testCovariate,
          extend = (end(sigRegions) - start(sigRegions) + 1)*2,
          addRegions = sigRegions,
          annoTrack = getAnnot(genome),
          lwd = 2,
          qval = FALSE,
          stat = FALSE,
          horizLegend = TRUE)
dev.off()

# Smoothed global and chromosomal methylation statistics  -----------------

dir.create("Global")

bs.filtered.bsseq %>%
  globalStats(testCovar = testCovariate,
              adjustCovar = adjustCovariate,
              matchCovar = matchCovariate) %>%
  openxlsx::write.xlsx("Global/smoothed_globalStats.xlsx") 

# PCAs of 20kb windows and CpG islands ------------------------------------

bs.filtered.bsseq %>%
  windowsPCA(goi = goi,
             group = bs.filtered.bsseq %>%
               pData() %>%
               dplyr::as_tibble() %>%
               dplyr::pull(!!testCovariate)
             ) %>% 
  ggplot2::ggsave("Global/Smoothed 20 Kb CpG Windows with CpG Islands.pdf",
                  plot = .,
                  device = NULL,
                  width = 11,
                  height = 8.5)

if(genome == "hg38" | genome == "hg19" | genome == "mm10" | genome == "mm9" | genome == "rn6"){
  bs.filtered.bsseq %>%
    CGiPCA(genome = genome, 
           group = bs.filtered.bsseq %>%
             pData() %>%
             dplyr::as_tibble() %>%
             dplyr::pull(!!testCovariate)
           ) %>% 
    ggplot2::ggsave("Global/Smoothed CpG Island Windows.pdf",
                    plot = .,
                    device = NULL,
                    width = 11,
                    height = 8.5)
}

bs.filtered.bsseq %>%
  densityPlot(goi = goi,
              group = bs.filtered.bsseq %>%
                pData() %>%
                dplyr::as_tibble() %>%
                dplyr::pull(!!testCovariate)
  ) %>% 
  ggplot2::ggsave("Global/Smoothed 20 Kb CpG Windows with CpG Islands Density Plot.pdf",
                  plot = .,
                  device = NULL,
                  width = 11,
                  height = 4)

# Heatmap -----------------------------------------------------------------

sigRegions %>%
  smoothPheatmap(bsseq = bs.filtered.bsseq,
                 testCovariate = testCovariate)

# CpG and genic annotations -----------------------------------------------

if(genome == "hg38" | genome == "hg19" | genome == "mm10" | genome == "mm9" | genome == "rn6"){
  annotateCpGs(sigRegions = sigRegions,
               regions = regions,
               genome = genome,
               saveAnnotations = T) %>%
    ggplot2::ggsave("DMRs/CpG_annotations.pdf",
                    plot = .,
                    device = NULL,
                    width = 8.5,
                    height = 11)

  annotateGenic(sigRegions = sigRegions,
                regions = regions,
                genome = genome,
                saveAnnotations = T) %>%
    ggplot2::ggsave("DMRs/generegion_annotations.pdf",
                    plot = .,
                    device = NULL,
                    width = 8.5,
                    height = 11)
}


# Gene symbol annotations -------------------------------------------------

cat("\n[DMRichR] Annotating DMRs with gene symbols \t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")

sigRegions %>%
  annotateRegions(TxDb = TxDb,
                  annoDb = annoDb) %T>%
  DMReport(regions = regions,
           bsseq = bs.filtered.bsseq,
           coverage = coverage,
           name = "DMReport") %>% 
  openxlsx::write.xlsx(file = "DMRs/DMRs_annotated.xlsx")

glue::glue("Annotating background regions with gene symbols...")
regions %>%
  annotateRegions(TxDb = TxDb,
                  annoDb = annoDb)
  openxlsx::write.xlsx(file = "DMRs/background_annotated.xlsx")
  
# Manhattan and Q-Q plots -------------------------------------------------

regionsAnno %>%
  manQQ()

# Gene ontology and pathway analyses  -------------------------------------

cat("\n[DMRichR] Performing gene ontology and pathway analyses \t", format(Sys.time(), "%d-%m-%Y %X"), "\n")

dir.create("Ontologies")

glue::glue("Running enrichR")
library(enrichR) # Needed or else "EnrichR website not responding"
#dbs <- listEnrichrDbs()
sigRegions %>%
  annotateRegions(TxDb = TxDb,
                  annoDb = annoDb) %>%  
  dplyr::select(geneSymbol) %>%
  purrr::flatten() %>%
  enrichR::enrichr(c("GO_Biological_Process_2018",
                     "GO_Cellular_Component_2018",
                     "GO_Molecular_Function_2018",
                     "KEGG_2019_Human",
                     "Panther_2016",
                     "Reactome_2016",
                     "RNA-Seq_Disease_Gene_and_Drug_Signatures_from_GEO")
          ) %T>%
  openxlsx::write.xlsx(file = "Ontologies/enrichr.xlsx") %>%
  GOplot(tool = "enrichR") %>%
  ggplot2::ggsave("Ontologies/enrichr_plot.pdf",
                  plot = .,
                  device = NULL,
                  height = 8.5,
                  width = 12)

glue::glue("Running rGREAT")
if(genome == "hg38" | genome == "hg19" | genome == "mm10" | genome == "mm9"){
  GREATjob <- sigRegions %>% 
    rGREAT::submitGreatJob(bg = regions,
                           species = genome,
                           request_interval = 1,
                           version = "4.0.4")
  
  glue::glue("Saving and plotting GREAT results...")
  GREATjob %>%
    rGREAT::getEnrichmentTables(category = "GO") %T>%
    openxlsx::write.xlsx(file = "Ontologies/GREAT_results.xlsx") %>% 
    GOplot(tool = "rGREAT") %>%
    ggplot2::ggsave("Ontologies/GREAT_plot.pdf",
           plot = .,
           device = NULL,
           height = 8.5,
           width = 12)
  
  pdf("Ontologies/GREAT_gene_associations_graph.pdf",
      height = 8.5,
      width = 11)
  par(mfrow = c(1, 3))
  res <- rGREAT::plotRegionGeneAssociationGraphs(GREATjob)
  dev.off()
  write.csv(as.data.frame(res),
            file = "Ontologies/GREATannotations.csv",
            row.names = F)
}

glue::glue("Running GOfuncR")
GOfuncR(sigRegions = sigRegions,
        regions = regions,
        genome = genome,
        n_randsets = 1000,
        upstream = 5000,
        downstream = 1000,
        annoDb = annoDb,
        TxDb = TxDb) %T>%
  openxlsx::write.xlsx("Ontologies/GOfuncR.xlsx") %>% 
  GOplot(tool = "GOfuncR") %>% 
  ggplot2::ggsave("Ontologies/GOfuncR_plot.pdf",
                  plot = .,
                  device = NULL,
                  height = 8.5,
                  width = 12)


# Machine learning --------------------------------------------------------

methylLearnOutput <- methylLearn(bsseq = bs.filtered.bsseq,
                                 regions = sigRegions,
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
GO_env <- ls(all = TRUE)[!(ls(all = TRUE) %in% bismark_env) &
                           !(ls(all = TRUE) %in% blocks_env) &
                           !(ls(all = TRUE) %in% DMRs_env) &
                           !(ls(all = TRUE) %in% bsseq_env)]
save(list = GO_env, file = "RData/GO.RData")
#load("RData/GO.RData")
                          
# End ---------------------------------------------------------------------

cat("\n[DMRichR] Summary \t\t\t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")

glue::glue("{length(sigRegions)} Significant DMRs \\
           ({round(sum(sigRegions$stat > 0) / length(sigRegions), digits = 2)*100}% hypermethylated, \\
           {round(sum(sigRegions$stat < 0) / length(sigRegions), digits = 2)*100}% hypomethylated) \\
           in {length(regions)} background regions \\
           from {nrow(bs.filtered)} CpGs assayed at {coverage}x coverage")

if(sum(blocks$pval < 0.05) > 0 & length(blocks) != 0){
glue::glue("{length(sigBlocks)} significant blocks of differential methylation in {length(blocks)} background blocks")
}

if(length(grep("genomecenter.ucdavis.edu", .libPaths())) == 0){sessionInfo()}
if(file.exists("Rplots.pdf")){file.remove("Rplots.pdf")}
rm(list = ls())
glue::glue("Done...")
quit(save = "no", status = 0, runLast = FALSE)
