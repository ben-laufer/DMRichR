#!/usr/bin/env Rscript

# DMRichR
# Ben Laufer

# R settings --------------------------------------------------------------

rm(list=ls())
options(scipen=999)
sink("DMRichR_log.txt", append = FALSE, split = TRUE)
#.libPaths("/share/lasallelab/programs/DMRichR/R_3.5")

# Functions ---------------------------------------------------------------

#' packageManage
#' @description Install package management
#' @export packageManage
packageManage <- function(){
  CRAN <- c("BiocManager", "remotes")
  new.packages <- CRAN[!(CRAN %in% installed.packages()[,"Package"])]
  if(length(new.packages)>0){
    install.packages(new.packages, repos ="https://cloud.r-project.org", quiet = TRUE)
  }
  stopifnot(suppressMessages(sapply(CRAN, require, character.only = TRUE)))
}

#' packageLoad
#' @description Install and load desired packages
#' @param packages Character string or desired packages
#' @export packageLoad
packageLoad <- function(packages = packages){
  new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)>0){
    library("BiocManager")
    new.packages <- gsub("ggbiplot", "vqv/ggbiplot", packages)
    new.packages <- gsub("DMRichR", "ben-laufer/DMRichR", packages)
    BiocManager::install(new.packages, ask = FALSE, quiet = TRUE)
  }
  stopifnot(suppressMessages(sapply(packages, require, character.only = TRUE)))
}

# Install and update ------------------------------------------------------

cat("\n[DMRichR] Installing and updating packages \t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
packageManage()
packageLoad(c("tidyverse", "dmrseq", "annotatr", "rGREAT", "enrichR", "ChIPseeker", "BiocParallel", "ggbiplot",
              "liftOver", "openxlsx", "CMplot", "optparse", "gplots", "RColorBrewer", "broom", "lsmeans", "DMRichR"))
suppressWarnings(BiocManager::valid(fix = TRUE, update = TRUE, ask = FALSE))

BiocManager::install("ben-laufer/DMRichR")
library(DMRichR)

# Global variables --------------------------------------------------------

cat("\n[DMRichR] Processing arguments from script \t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")

option_list <- list(
  make_option(c("-g", "--genome"), type = "character", default = NULL,
              help = "Choose a genome (hg38, mm10, rn6, rheMac8) [required]"),
  make_option(c("-x", "--coverage"), type = "integer", default = 1,
              help = "Choose a CpG coverage cutoff [default = %default]"),
  make_option(c("-t", "--testCovariate"), type = "character", default = NULL,
              help = "Choose a test covariate [required]"),
  make_option(c("-a", "--adjustCovariate"), type = "character", default = NULL,
              help = "Choose covariates to directly adjust [default = NULL]"),
  make_option(c("-m", "--matchCovariate"), type = "character", default = NULL,
              help = "Choose covariate to balance permutations [default = NULL]"),
  make_option(c("-c", "--cores"), type = "integer", default = 8,
              help = "Choose number of cores [default = %default]")
)
opt <- parse_args(OptionParser(option_list=option_list))

message("Assigning arguments to global variables...")
# Check for requirements
stopifnot(!is.null(opt$genome))
stopifnot(!is.null(opt$testCovariate))
# Assign
genome <- as.character(opt$genome)
coverage <- as.numeric(opt$coverage)
testCovariate <- as.character(opt$testCovariate)
if(!is.null(opt$adjustCovariate)){
  adjustCovariate <- opt$adjustCovariate %>% strsplit(";") %>% unlist() %>% as.character()
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
cat(paste("genome =", genome), "\n")
cat(paste("coverage =", coverage), "\n")
cat(paste("testCovariate =", testCovariate), "\n")
cat(paste("adjustCovariate =", adjustCovariate), "\n")
cat(paste("matchCovariate =", matchCovariate), "\n")
cat(paste("cores =", cores), "\n")

# Setup Annotation Databases ----------------------------------------------

cat("\n[DMRichR] Selecting annotation databases \t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")

if(genome == "hg38"){
  packages <- c("BSgenome.Hsapiens.UCSC.hg38", "TxDb.Hsapiens.UCSC.hg38.knownGene", "org.Hs.eg.db")
}else if(genome == "mm10"){
  packages <- c("BSgenome.Mmusculus.UCSC.mm10", "TxDb.Mmusculus.UCSC.mm10.knownGene", "org.Mm.eg.db")
}else if(genome == "rheMac8"){
  packages <- c("BSgenome.Mmulatta.UCSC.rheMac8", "TxDb.Mmulatta.UCSC.rheMac8.refGene", "org.Mmu.eg.db")
}else if(genome == "rn6"){
  packages <- c("BSgenome.Rnorvegicus.UCSC.rn6", "TxDb.Rnorvegicus.UCSC.rn6.refGene", "org.Rn.eg.db")
}else{
  stop(paste(genome, "is not suppourted, please choose either hg38, mm10, rheMac8, or rn6 [Case Sensitive]"))
}

packageLoad(packages)

if(genome == "hg38"){
  goi <- BSgenome.Hsapiens.UCSC.hg38; TxDb <- TxDb.Hsapiens.UCSC.hg38.knownGene; annoDb <- "org.Hs.eg.db"
}else if(genome == "mm10"){
  goi <- BSgenome.Mmusculus.UCSC.mm10; TxDb <- TxDb.Mmusculus.UCSC.mm10.knownGene; annoDb <- "org.Mm.eg.db"
}else if(genome == "rheMac8"){
  goi <- BSgenome.Mmulatta.UCSC.rheMac8; TxDb <- TxDb.Mmulatta.UCSC.rheMac8.refGene; annoDb <- "org.Mmu.eg.db"
}else if(genome == "rn6"){
  goi <- BSgenome.Rnorvegicus.UCSC.rn6; TxDb <- TxDb.Rnorvegicus.UCSC.rn6.refGene; annoDb <- "org.Rn.eg.db"
}else{
  stop(paste(genome, "is not suppourted, please choose either hg38, mm10, rheMac8, or rn6 [Case Sensitive]"))
}

# Load and process samples ------------------------------------------------

name <- gsub( "_.*$","", list.files(path=getwd(), pattern="*.txt.gz"))

bs.filtered <- processBismark(files = list.files(path=getwd(), pattern="*.txt.gz"),
                              meta = read.csv("sample_info.csv", header = TRUE),
                              groups = testCovariate,
                              Cov = coverage,
                              mc.cores = cores)

message("Saving Rdata...")
bismark_env <- ls(all = TRUE)
save(list = bismark_env, file = "bismark.RData")
#load("bismark.RData")

# Distribtuion plots ------------------------------------------------------

#cat("\n[DMRichR] Plotting Empirical Distribution of CpGs \t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
#pdf("Filtered_CpG_Methylation_Distributions.pdf", height = 7.50, width = 11.50)
#plotEmpiricalDistribution(bs.filtered, testCovariate = testCovariate)
#plotEmpiricalDistribution(bs.filtered, testCovariate = testCovariate, type = "Cov", bySample = TRUE)
#dev.off()

# DMRs --------------------------------------------------------------------

cat("\n[DMRichR] Testing for DMRs with dmrseq \t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
# Reproducible permutations (change and record seed for different datasets to avoid any potential random bias)
set.seed(5)
#.Random.seed
# Only use 1 core, multicore is slower due to forking problem (for now?)
register(MulticoreParam(1))
regions <- dmrseq(bs=bs.filtered,
                  cutoff = 0.05,
                  minNumRegion = 5,
                  maxPerms = 10,
                  testCovariate = testCovariate,
                  adjustCovariate = adjustCovariate,
                  matchCovariate = matchCovariate)

message("Selecting significant DMRs...")
if(sum(regions$qval < 0.05) < 100 & sum(regions$pval < 0.05) != 0){
  sigRegions <- regions[regions$pval < 0.05,]
}else if(sum(regions$qval < 0.05) >= 100){
  sigRegions <- regions[regions$qval < 0.05,]
}else if(sum(regions$pval < 0.05) == 0){
  stop("No significant DMRs detected")
  }
cat(paste(round(sum(sigRegions$stat > 0) / length(sigRegions), digits = 2)*100, "% of DMRs are hypermethylated", sep =""), "\n")

message("Plotting DMR pie chart...")
pie <- (table(sigRegions$stat < 0))
names(pie) <- c("Hypermethylated", "Hypomethylated")
pdf("HypervsHypo_pie.pdf", height = 8.5, width = 11)
pie(pie,
    labels = c(paste(pie[1], "Hypermethylated", sep = " "), paste(pie[2], "Hypomethylated", sep = " ")),
    col = c("Red", "Blue"))
dev.off()

#message("Extracing raw differnces for DMRs...")
# Does not work anymore, is it from new bsseq:read.bismark() changes to index? Error: subscript contains out-of-bounds ranges
#rawDiff <- meanDiff(bs.filtered,
#                    dmrs = regions,
#                    testCovariate = testCovariate)
#sigRawDiff <- meanDiff(bs.filtered,
#                       dmrs = sigRegions,
#                       testCovariate = testCovariate)
#regions$RawDiff <- rawDiff
#sigRegions$RawDiff <- sigRawDiff

message("Calculating average percent differences...") 
regions$percentDifference <- round(regions$beta/pi * 100)
sigRegions$percentDifference <- round(sigRegions$beta/pi *100)

message("Exporting DMR and background region information...")
gr2csv(regions, "backgroundRegions.csv")
gr2csv(sigRegions, "DMRs.csv")
gr2bed(regions, "backgroundRegions.bed")
gr2bed(sigRegions, "DMRs.bed")

message("Annotating and plotting...")
pdf("DMRs.pdf", height = 7.50, width = 11.50)
annoTrack <- getAnnot(genome)
plotDMRs(bs.filtered,
         regions = sigRegions,
         testCovariate = testCovariate,
         annoTrack = annoTrack,
         qval = F)
dev.off()

message("Saving Rdata...")
DMRs_env <- ls(all = TRUE)[!(ls(all = TRUE) %in% bismark_env)]
save(list = DMRs_env, file = "DMRs.RData")
#load("DMRs.RData")

# Individual smoothed values ----------------------------------------------

cat("\n[DMRichR] Smoothing individual methylation values \t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
bs.filtered.bsseq <- BSmooth(bs.filtered,
                             BPPARAM = MulticoreParam(workers = cores, progressbar = TRUE))
bs.filtered.bsseq

message("Extracting values for WGCNA...")
indiv_smoothed_table <- getSmooth(bsseq = bs.filtered.bsseq,
                                  regions = regions,
                                  out = "background_region_individual_smoothed_methylation.txt")

message("Saving Rdata...")
bsseq_env <- ls(all = TRUE)[!(ls(all = TRUE) %in% bismark_env) &
                              !(ls(all = TRUE) %in% DMRs_env)]
save(list = bsseq_env, file = "bsseq.RData")
#load("bsseq.RData")

# Global methylation ------------------------------------------------------

bs.filtered.bsseq %>%
  getGlobal() %>%
  smoothANOVA() %>%
  write.xlsx("smoothed_global_methylation_stats.xlsx")

# Chromosomal methylation -------------------------------------------------

bs.filtered.bsseq %>%
  getChrom() %>%
  smoothANOVA() %>%
  write.xlsx("smoothed_global_chromosomal_methylation_stats.xlsx")

# PCA of 20 kb windows with CGi -------------------------------------------

cat("\n[DMRichR] 20 kb windows \t\t\t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
message("Creating 20 kb windows")
chrSizes <- seqlengths(goi)
windows <- tileGenome(chrSizes,
                      tilewidth = 2e4,
                      cut.last.tile.in.chrom = TRUE)
windows <- GenomeInfoDb::keepStandardChromosomes(windows, pruning.mode = "coarse")
windows

message("Extracting values for 20 kb windows...")
windows_smoothed_table <- getSmooth(bsseq = bs.filtered.bsseq,
                                    regions = windows,
                                    out = "20kb_smoothed_windows.txt")

message("Tidying 20 kb window data...")
meth_reorder <- na.omit(windows_smoothed_table[,c(6:length(windows_smoothed_table))])
data <- t(as.matrix(meth_reorder))
# Fix for columns of no variance
#data2 <- as.data.frame(data)
#data3 <- data2[,apply(data2, 2, var, na.rm=TRUE) != 0]
# Titles
stopifnot(sampleNames(bs.filtered.bsseq) == colnames(meth_reorder))
group <- bs.filtered.bsseq %>% pData() %>% as.tibble() %>% pull(!!testCovariate)

message("20 kb window PCA...")
PCA(data, "Smoothed 20 Kb CpG Windows with CpG Islands")

# PCA of CGi windows ------------------------------------------------------

if(genome == "hg38" | genome == "mm10" | genome == "rn6"){
  cat("\n[DMRichR] CGi windows \t\t\t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
  message("Obtaining CGi annotations...")
  CGi <- build_annotations(genome = genome, annotations = paste(genome,"_cpg_islands", sep = ""))
  CGi <- GenomeInfoDb::keepStandardChromosomes(CGi, pruning.mode = "coarse")
  CGi

  message("Extracting values for CGi windows...")
  CGi_smoothed_table <- getSmooth(bsseq = bs.filtered.bsseq,
                                  regions = CGi,
                                  out = "CGi_smoothed_windows.txt")

  message("Tidying CGi window data...")
  meth_reorder <- na.omit(CGi_smoothed_table[,c(11:length(CGi_smoothed_table))])
  data <- t(as.matrix(meth_reorder))
  # Fix for columns of no variance
  #data2 <- as.data.frame(data)
  #data3 <- data2[,apply(data2, 2, var, na.rm=TRUE) != 0]
  # Titles
  stopifnot(sampleNames(bs.filtered.bsseq) == colnames(meth_reorder))
  group <- bs.filtered.bsseq %>% pData() %>% as.tibble() %>% pull(!!testCovariate)

  message("CGi window PCA...")
  PCA(data, "Smoothed CpG Island Windows")
}

# Heatmap -----------------------------------------------------------------

smoothHeatmap(regions = sigRegions,
              bsseq = bs.filtered.bsseq,
              groups = bs.filtered.bsseq %>% pData() %>% as.tibble() %>% pull(!!testCovariate),
              out = "sig_individual_smoothed_DMR_methylation.txt")

# Prepare files for enrichment analyses -----------------------------------

if(genome == "hg38" | genome == "mm10" | genome == "rn6"){
  cat("\n[DMRichR] Preparing DMRs files for annotations \t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
  external <- sigRegions
  for (i in 1:length(external)){
    if(external$stat[i] > 0){
      external$direction[i] <- "Hypermethylated"
    }else if(external$stat[i] < 0){
      external$direction[i] <- "Hypomethylated"
    }else{
      stop("Annotation problem")
    }}
  externalOut <- as.data.frame(external)
  dir.create("GAT")
  df2bed(externalOut[,c(1:3,16)], "GAT/DMRs.bed")

  message("Preparing DMRs for HOMER...")
  dir.create("HOMER")
  gr2bed((external[,c(1:3)])[external$direction == "Hypermethylated",], "HOMER/DMRs_hyper.bed")
  gr2bed((external[,c(1:3)])[external$direction == "Hypomethylated",], "HOMER/DMRs_hypo.bed")

  message("Preparing background regions for annotations...")
  external_bg <- regions
  for (i in 1:length(external_bg)){
    if(external_bg$stat[i] > 0){
      external_bg$direction[i] <- "Hypermethylated"
    }else if(external_bg$stat[i] < 0){
      external_bg$direction[i] <- "Hypomethylated"
    }else{
      stop("Annotation problem")
    }}

  # CpG Annotations ---------------------------------------------------------

  cat("\n[DMRichR] Building CpG annotations \t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
  annotations <- build_annotations(genome = genome, annotations = paste(genome,"_cpgs", sep=""))
  annotations <- GenomeInfoDb::keepStandardChromosomes(annotations, pruning.mode = "coarse")

  message("Annotating DMRs...")
  dm_annotated_CpG <- annotate_regions(
    regions = external,
    annotations = annotations,
    ignore.strand = TRUE,
    quiet = FALSE)

  message("Annotating background regions...")
  background_annotated_CpG <- annotate_regions(
    regions = external_bg,
    annotations = annotations,
    ignore.strand = TRUE,
    quiet = FALSE)

  message("Saving files for GAT...")
  CpGs <- as.data.frame(annotations)
  CpGs <- CpGs[!grepl("_", CpGs$seqnames), ]
  table(CpGs$seqnames)
  df2bed(CpGs[, c(1:3,10)], paste("GAT/", genome, "CpG.bed", sep = ""))

  message("Preparing CpG annotation plot...")
  x_order <- c('Hypermethylated','Hypomethylated')
  fill_order <- c(
    paste(genome,"_cpg_islands",sep=""),
    paste(genome,"_cpg_shores",sep=""),
    paste(genome,"_cpg_shelves",sep=""),
    paste(genome,"_cpg_inter",sep=""))

  CpG_bar <- plot_categorical(
    annotated_regions = dm_annotated_CpG,
    annotated_random = background_annotated_CpG,
    x = 'direction',
    fill = 'annot.type',
    x_order = x_order,
    fill_order = fill_order,
    position='fill',
    plot_title = '',
    legend_title = 'Annotations',
    x_label = '',
    y_label = 'Proportion') +
    scale_x_discrete(labels=c("All", "Hypermethylated", "Hypomethylated", "Background")) +
    scale_y_continuous(expand=c(0,0)) +
    theme_classic() +
    theme(axis.text = element_text(size = 25),
          axis.title = element_text(size = 25),
          strip.text = element_text(size = 25),
          legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave("CpG_annotations.pdf", plot = CpG_bar, device = NULL, width = 8.5, height = 11)

  # Gene Annotations --------------------------------------------------------

  cat("\n[DMRichR] Building gene region annotations \t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
  annotations <- build_annotations(genome = genome, annotations = c(paste(genome,"_basicgenes", sep = ""),
                                                                    paste(genome,"_genes_intergenic", sep = ""),
                                                                    paste(genome,"_genes_intronexonboundaries", sep = ""),
                                                                    if(genome == "hg38" | genome == "mm10"){paste(genome,"_enhancers_fantom", sep = "")}))
  annotations <- GenomeInfoDb::keepStandardChromosomes(annotations, pruning.mode = "coarse")

  message("Saving files for GAT...")
  annoFile <- as.data.frame(annotations)
  annoFile <- annoFile[!grepl("_", annoFile$seqnames) ,]
  table(annoFile$seqnames)
  annoFile <- annoFile[, c(1:3,10)]

  if(genome == "hg38" | genome == "mm10"){
    gr2bed(annoFile[annoFile$type == paste(genome,"_enhancers_fantom", sep = ""), ], "GAT/enhancers.bed")}
  gr2bed(annoFile[annoFile$type == paste(genome,"_genes_promoters", sep = ""), ], "GAT/promoters.bed")
  gr2bed(annoFile[annoFile$type == paste(genome,"_genes_introns", sep = ""), ], "GAT/introns.bed")
  gr2bed(annoFile[annoFile$type == paste(genome,"_genes_intronexonboundaries", sep = ""), ], "GAT/boundaries.bed")
  gr2bed(annoFile[annoFile$type == paste(genome,"_genes_intergenic", sep = ""), ], "GAT/intergenic.bed")
  gr2bed(annoFile[annoFile$type == paste(genome,"_genes_exons", sep = ""), ], "GAT/exons.bed")
  gr2bed(annoFile[annoFile$type == paste(genome,"_genes_5UTRs", sep = ""), ], "GAT/fiveUTRs.bed")
  gr2bed(annoFile[annoFile$type == paste(genome,"_genes_3UTRs", sep = ""), ], "GAT/threeUTRs.bed")
  gr2bed(annoFile[annoFile$type == paste(genome,"_genes_1to5kb", sep = ""), ], "GAT/onetofivekb.bed")

  message("Annotating DMRs...")
  dm_annotated <- annotate_regions(
    regions = external,
    annotations = annotations,
    ignore.strand = TRUE,
    quiet = FALSE)

  message("Annotating background regions...")
  background_annotated <- annotate_regions(
    regions = external_bg,
    annotations = annotations,
    ignore.strand = TRUE,
    quiet = FALSE)

  message("Preparing CpG annotation plot...")
  x_order <- c('Hypermethylated','Hypomethylated')
  fill_order <- c(
    if(genome == "hg38" | genome == "mm10"){paste(genome, "_enhancers_fantom", sep = "")},
    paste(genome,"_genes_1to5kb", sep = ""),
    paste(genome,"_genes_promoters", sep = ""),
    paste(genome,"_genes_5UTRs", sep = ""),
    paste(genome,"_genes_exons", sep = ""),
    paste(genome,"_genes_intronexonboundaries", sep = ""),
    paste(genome,"_genes_introns", sep = ""),
    paste(genome,"_genes_3UTRs", sep = ""),
    paste(genome,"_genes_intergenic", sep = ""))

  gene_bar <- plot_categorical(
    annotated_regions = dm_annotated,
    annotated_random = background_annotated,
    x = 'direction',
    fill = 'annot.type',
    x_order = x_order,
    fill_order = fill_order,
    position = 'fill',
    plot_title = '',
    legend_title = 'Annotations',
    x_label = '',
    y_label = 'Proportion') +
    scale_x_discrete(labels=c("All", "Hypermethylated", "Hypomethylated", "Background")) +
    scale_y_continuous(expand=c(0,0)) +
    theme_classic() +
    theme(axis.text = element_text(size = 25),
          axis.title = element_text(size = 25),
          strip.text = element_text(size = 25),
          legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave("generegion_annotations.pdf", plot = gene_bar, device = NULL, width = 8.5, height = 11)
}

# GREAT -------------------------------------------------------------------

if(genome=="hg38"){
  cat("\n[DMRichR] Obtaining liftOver information for GREAT \t", format(Sys.time(), "%d-%m-%Y %X"))
  path <- system.file(package="liftOver", "extdata", "hg38ToHg19.over.chain")
  ch <- import.chain(path)
  ch

  message("liftOver DMRs...")
  seqlevelsStyle(sigRegions) <- "UCSC"
  sigRegions_liftOver <- liftOver(sigRegions, ch)
  class(sigRegions_liftOver)
  sigRegions_liftOver <- unlist(sigRegions_liftOver)
  length(sigRegions) - length(sigRegions_liftOver)

  message("liftOver background regions...")
  seqlevelsStyle(regions) <- "UCSC"
  regions_liftOver <- liftOver(regions, ch)
  class(regions_liftOver)
  regions_liftOver <- unlist(regions_liftOver)
  length(regions) - length(regions_liftOver)
}

if(genome == "hg38" | genome == "mm10"){
  cat("\n[DMRichR] Submitting to GREAT\ \t\t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
  if(genome == "hg38"){
    gr <- sigRegions_liftOver; bg <- regions_liftOver; species <- "hg19"
  }else if(genome == "mm10"){
    gr <- sigRegions; bg <- regions; species <- "mm10"
  }else{
    stop(paste(genome, "is not suppourted for GREAT, please choose either hg38 or mm10 [Case Sensitive]"))
  }

  job <- submitGreatJob(gr,
                        bg = bg,
                        species = species)
  job
  #availableCategories(job)
  #availableOntologies(job)
  tb <- getEnrichmentTables(job, category = c("GO", "Pathway Data"))

  message("Saving GREAT enrichment results...")
  write.xlsx(tb, file = "GREAT_results.xlsx", sep="")

  message("Plotting GREAT results...")
  pdf("GREAT.pdf", height = 8.5, width = 11)
  par(mfrow = c(1, 3))
  res <- plotRegionGeneAssociationGraphs(job)
  dev.off()

  message("Saving GREAT annotations...")
  res
  write.csv(as.data.frame(res), file="GREATannotations.csv", row.names = F)
}

# Chipseeker --------------------------------------------------------------

cat("\n[DMRichR] Annotating DMRs with gene symbols \t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
peakAnno <- annotatePeak(sigRegions,
                         TxDb = TxDb,
                         annoDb = annoDb,
                         overlap = "all")

message("Annotating background regions with gene symbols...")
backgroundAnno <- annotatePeak(regions,
                               TxDb = TxDb,
                               annoDb = annoDb,
                               overlap = "all")

message("Upset Plot of genic features...")
pdf("Upset.pdf", onefile = FALSE, height = 8.5, width = 11)
upsetplot(peakAnno, vennpie = TRUE)
dev.off()

message("Saving gene annotations...")
annotations <- as.data.frame(peakAnno)
annotations[2] <- annotations[2] - 1
annotations <- annotations[-c(5,17:19)]
write.xlsx(annotations, file = "DMRs_annotated.xlsx", sep= "")

background_annotations <- as.data.frame(backgroundAnno)
background_annotations[2] <- background_annotations[2] - 1
background_annotations <- background_annotations[-c(5,17:19)]
write.xlsx(background_annotations, file = "background_annotated.xlsx", sep= "")

# CMplot ------------------------------------------------------------------

cat("\n[DMRichR] Manhattan and QQ plots \t\t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
message("Tidying for Manhattan and QQ plots", format(Sys.time(), "%d-%m-%Y %X"), "\n")
Manhattan <- as.data.frame(sort(as.GRanges(backgroundAnno), ignore.strand=TRUE))[c("SYMBOL","seqnames", "start", "pval")]
Manhattan$seqnames <- substring(Manhattan$seqnames, 4)
cols = gg_color_hue(2)

message("Generating Manhattan and QQ plots...")
CMplot(Manhattan,
       col = cols,
       plot.type = c("m","q"),
       LOG10 = TRUE,
       ylim = NULL,
       threshold = 0.05, #c(1e-6,1e-4),
       threshold.lty = c(1,2),
       threshold.lwd = c(1,1),
       threshold.col = c("black","grey"),
       cex = 0.5,
       cex.axis = 0.7,
       amplify = FALSE,
       chr.den.col = brewer.pal(9, "YlOrRd"),
       bin.size = 1e6,
       bin.max = 100,
       signal.col = c("red","green"),
       signal.cex = c(1,1),
       signal.pch = c(19,19),
       file = "pdf",
       memo = "")

# Enrichr -----------------------------------------------------------------

cat("\n[DMRichR] Running Enrichr \t\t\t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
# Check available databases
#dbs <- listEnrichrDbs()
dbs <- c("GO_Biological_Process_2018",
         "GO_Cellular_Component_2018",
         "GO_Molecular_Function_2018",
         "KEGG_2016",
         "Panther_2016",
         "Reactome_2016",
         "RNA-Seq_Disease_Gene_and_Drug_Signatures_from_GEO")
GO <- enrichr(annotations$SYMBOL, dbs)
write.xlsx(GO, file = "enrichr.xlsx", sep="")

message("Saving RData...")
GO_env <- ls(all = TRUE)[!(ls(all = TRUE) %in% bismark_env) &
                           !(ls(all = TRUE) %in% DMRs_env) &
                           !(ls(all = TRUE) %in% bsseq_env)]
save(list = GO_env, file = "GO.RData")
#load("GO.RData")

# Blocks ------------------------------------------------------------------

cat("\n[DMRichR] Testing for large blocks (PMDs/HMDs) \t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")

register(MulticoreParam(1))
blocks <- dmrseq(bs = bs.filtered,
                 cutoff = 0.05,
                 maxPerms = 10,
                 testCovariate=testCovariate,
                 adjustCovariate = adjustCovariate,
                 matchCovariate = matchCovariate,
                 block = TRUE,
                 minInSpan = 500,
                 bpSpan = 5e4,
                 maxGapSmooth = 1e6,
                 maxGap = 5e3)

message("Selecting significant blocks...")

if(sum(blocks$qval < 0.05) == 0 & sum(blocks$pval < 0.05) != 0){
  sigBlocks <- blocks[blocks$pval < 0.05,]
}else if(sum(blocks$qval < 0.05) >= 1){
  sigBlocks <- blocks[blocks$qval < 0.05,]
}else if(sum(blocks$pval < 0.05) == 0){
  print("No significant blocks detected")
}

message("Exporting block and background information...")
gr2csv(blocks, "backgroundBlocks.csv")
gr2bed(blocks, "backgroundBlocks.bed")
if(sum(blocks$pval < 0.05) > 0){
  gr2csv(sigBlocks, "blocks.csv")
  gr2bed(sigBlocks, "blocks.bed")
}

if(length(sigBlocks) > 0){
  message("Annotating and plotting blocks...")
  pdf("Blocks.pdf", height = 7.50, width = 11.50)
  annoTrack <- getAnnot(genome)
  plotDMRs(bs.filtered,
           regions = sigBlocks,
           testCovariate = testCovariate,
           annoTrack = annoTrack,
           qval = FALSE)
  dev.off()
}

message("Saving RData...")
blocks_env <- ls(all = TRUE)[!(ls(all = TRUE) %in% bismark_env) &
                               !(ls(all = TRUE) %in% DMRs_env) &
                               !(ls(all = TRUE) %in% bsseq_env) &
                               !(ls(all = TRUE) %in% GO_env)]
save(list = blocks_env, file = "Blocks.RData")
#load("Blocks.RData")

# End ---------------------------------------------------------------------

cat("\n[DMRichR] Finishing \t\t\t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")

sessionInfo()
rm(list = ls())
message("Done...")

