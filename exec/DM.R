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

if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")
devtools::install_github("ben-laufer/DMRichR")

DMRichR::packageLoad(c("tidyverse", "dmrseq", "annotatr", "rGREAT", "enrichR", "ChIPseeker", "BiocParallel", "ggbiplot",
                       "liftOver", "openxlsx", "CMplot", "optparse", "gplots", "RColorBrewer", "broom", "lsmeans", "glue",
                       "caret", "e1071", "randomForest", "randomForestExplainer", "gt", "DMRichR"))

# Global variables --------------------------------------------------------

cat("\n[DMRichR] Processing arguments from script \t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")

option_list <- list(
  make_option(c("-g", "--genome"), type = "character", default = NULL,
              help = "Choose a genome (hg38, mm10, rn6, rheMac8) [required]"),
  make_option(c("-x", "--coverage"), type = "integer", default = 1,
              help = "Choose a CpG coverage cutoff [default = %default]"),
  make_option(c("-s", "--perGroup"), type = "double", default = 1,
              help = "Choose the percent [values from 0 to 1] of samples in a 2 factor group for CpG coverage cutoff [default = %default]"),
  make_option(c("-n", "--minCpGs"), type = "integer", default = 5,
              help = "Choose the minimum number of CpGs for a DMR [default = %default]"),
  make_option(c("-p", "--maxPerms"), type = "integer", default = 10,
              help = "Choose the number of permutations for DMR and block analyses [default = %default]"),
  make_option(c("-o", "--cutoff"), type = "double", default = 0.05,
              help = "Choose the cutoff value [from 0 to 1] for the single CpG coefficient utilized to discover testable background regions [default = %default]"),
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

glue::glue("Assigning arguments to global variables...")
# Check for requirements
stopifnot(!is.null(opt$genome))
stopifnot(!is.null(opt$testCovariate))
# Assign
genome <- as.character(opt$genome)
coverage <- as.numeric(opt$coverage)
perGroup <- as.numeric(opt$perGroup)
minCpGs <- as.numeric(opt$minCpGs)
maxPerms <- as.numeric(opt$maxPerms)
cutoff <- as.numeric(opt$cutoff)
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
glue("genome = {genome}")
glue("coverage = {coverage}")
glue("perGroup = {perGroup}")
glue("minCpGs = {minCpGs}")
glue("maxPerms = {maxPerms}")
glue("cutoff = {cutoff}")
glue("testCovariate = {testCovariate}")
glue("adjustCovariate = {adjustCovariate}")
glue("matchCovariate = {matchCovariate}")
glue("cores = {cores}")

# Setup annotation databases ----------------------------------------------

cat("\n[DMRichR] Selecting annotation databases \t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")

packages <- dplyr::case_when(genome == "hg38" ~ c("BSgenome.Hsapiens.UCSC.hg38", "TxDb.Hsapiens.UCSC.hg38.knownGene", "org.Hs.eg.db"),
                             genome == "mm10" ~ c("BSgenome.Mmusculus.UCSC.mm10", "TxDb.Mmusculus.UCSC.mm10.knownGene", "org.Mm.eg.db"),
                             genome == "rheMac8" ~ c("BSgenome.Mmulatta.UCSC.rheMac8", "TxDb.Mmulatta.UCSC.rheMac8.refGene", "org.Mmu.eg.db"),
                             genome == "rn6" ~ c("BSgenome.Rnorvegicus.UCSC.rn6", "TxDb.Rnorvegicus.UCSC.rn6.refGene", "org.Rn.eg.db")
                             )

packageLoad(packages)

if(genome == "hg38"){
  goi <- BSgenome.Hsapiens.UCSC.hg38
  TxDb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  annoDb <- "org.Hs.eg.db"
}else if(genome == "mm10"){
  goi <- BSgenome.Mmusculus.UCSC.mm10
  TxDb <- TxDb.Mmusculus.UCSC.mm10.knownGene
  annoDb <- "org.Mm.eg.db"
}else if(genome == "rheMac8"){
  goi <- BSgenome.Mmulatta.UCSC.rheMac8
  TxDb <- TxDb.Mmulatta.UCSC.rheMac8.refGene
  annoDb <- "org.Mmu.eg.db"
}else if(genome == "rn6"){
  goi <- BSgenome.Rnorvegicus.UCSC.rn6
  TxDb <- TxDb.Rnorvegicus.UCSC.rn6.refGene
  annoDb <- "org.Rn.eg.db"
}else{
  stop(glue("{genome} is not supported, please choose either hg38, mm10, rheMac8, or rn6 [Case Sensitive]"))
}

# Load and process samples ------------------------------------------------

name <- gsub( "_.*$","", list.files(path = getwd(), pattern = "*.txt.gz"))

bs.filtered <- processBismark(files = list.files(path = getwd(), pattern = "*.txt.gz"),
                              meta = read.xlsx("sample_info.xlsx", colNames = TRUE) %>% mutate_if(is.character,as.factor),
                              groups = testCovariate,
                              Cov = coverage,
                              mc.cores = cores,
                              per.Group = perGroup)

glue::glue("Saving Rdata...")
bismark_env <- ls(all = TRUE)
save(list = bismark_env, file = "bismark.RData")
#load("bismark.RData")

# Distribution plots ------------------------------------------------------

#cat("\n[DMRichR] Plotting Empirical Distribution of CpGs \t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
#pdf("Filtered_CpG_Methylation_Distributions.pdf", height = 7.50, width = 11.50)
#plotEmpiricalDistribution(bs.filtered, testCovariate = testCovariate)
#plotEmpiricalDistribution(bs.filtered, testCovariate = testCovariate, type = "Cov", bySample = TRUE)
#dev.off()

# Background --------------------------------------------------------------

cat("\n[DMRichR] Getting bsseq background regions \t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
background <- getBackground(bs.filtered,
                            minNumRegion = minCpGs,
                            maxGap = 1000)
write.table(background,
            file = "bsseq_background.csv",
            sep = ",",
            quote = FALSE,
            row.names = FALSE)

# DMRs --------------------------------------------------------------------

cat("\n[DMRichR] Testing for DMRs with dmrseq \t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
start_time <- Sys.time()

# Reproducible permutations (change and record seed for different datasets to avoid any potential random bias)
set.seed(5)
#.Random.seed

# More cores increases smoothing time but decreases scoring time, so this is my attempt at balancing it
glue::glue("Determining parallelization...") 
if(cores >= 4){
  BPPARAM <- BiocParallel::MulticoreParam(workers = ceiling(cores/4))
  glue::glue("Parallel processing will be used with {ceiling(cores/4)} cores")
}else if(cores < 4){
  BPPARAM <- BiocParallel::MulticoreParam(workers = 1)
  glue::glue("Parallel processing will not be used")
}
register(BPPARAM)

regions <- dmrseq(bs = bs.filtered,
                  cutoff = cutoff,
                  minNumRegion = minCpGs,
                  maxPerms = maxPerms,
                  testCovariate = testCovariate,
                  adjustCovariate = adjustCovariate,
                  matchCovariate = matchCovariate)

glue::glue("Selecting significant DMRs...", "\n")
if(sum(regions$qval < 0.05) < 100 & sum(regions$pval < 0.05) != 0){
  sigRegions <- regions[regions$pval < 0.05,]
}else if(sum(regions$qval < 0.05) >= 100){
  sigRegions <- regions[regions$qval < 0.05,]
}else if(sum(regions$pval < 0.05) == 0){
  stop(glue("No significant DMRs detected in {length(regions)} background regions"))
  }

if(sum(sigRegions$stat > 0) > 0 & sum(sigRegions$stat < 0) > 0){
  glue::glue("{length(sigRegions)} Significant DMRs \\
             ({round(sum(sigRegions$stat > 0) / length(sigRegions), digits = 2)*100}% hypermethylated, \\
             {round(sum(sigRegions$stat < 0) / length(sigRegions), digits = 2)*100}% hypomethylated) \\
             in {length(regions)} background regions \\
             from {nrow(bs.filtered)} CpGs assayed at {coverage}x coverage")
}

glue::glue("Calculating average percent differences...") 
regions$percentDifference <- round(regions$beta/pi * 100)
sigRegions$percentDifference <- round(sigRegions$beta/pi *100)

glue::glue("Exporting DMR and background region information...")
gr2csv(regions, "backgroundRegions.csv")
gr2csv(sigRegions, "DMRs.csv")
gr2bed(regions, "backgroundRegions.bed")
gr2bed(sigRegions, "DMRs.bed")

# glue::glue("Annotating and plotting...")
# pdf("DMRs.pdf", height = 7.50, width = 11.50)
# annoTrack <- getAnnot(genome)
# plotDMRs(bs.filtered,
#          regions = sigRegions,
#          testCovariate = testCovariate,
#          annoTrack = annoTrack,
#          qval = F)
# dev.off()

glue::glue("Saving Rdata...")
DMRs_env <- ls(all = TRUE)[!(ls(all = TRUE) %in% bismark_env)]
save(list = DMRs_env, file = "DMRs.RData")
#load("DMRs.RData")

glue::glue("DMR timing...")
end_time <- Sys.time()
end_time - start_time

# Individual smoothed values ----------------------------------------------

cat("\n[DMRichR] Smoothing individual methylation values \t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
start_time <- Sys.time()

bs.filtered.bsseq <- BSmooth(bs.filtered,
                             BPPARAM = MulticoreParam(workers = ceiling(cores/3), progressbar = FALSE))

# Drop chrY in Rat only due to poor quality (some CpGs in females map to Y)
if(genome == "rn6"){
  bs.filtered.bsseq <- dropSeqlevels(bs.filtered.bsseq, "chrY", pruning.mode = "coarse")
  seqlevels(bs.filtered.bsseq)
}

bs.filtered.bsseq

glue::glue("Extracting values for WGCNA...")
indiv_smoothed_table <- getSmooth(bsseq = bs.filtered.bsseq,
                                  regions = regions,
                                  out = "background_region_individual_smoothed_methylation.txt")

glue::glue("Saving Rdata...")
bsseq_env <- ls(all = TRUE)[!(ls(all = TRUE) %in% bismark_env) &
                              !(ls(all = TRUE) %in% DMRs_env)]
save(list = bsseq_env, file = "bsseq.RData")
#load("bsseq.RData")

glue::glue("Individual smoothing timing...")
end_time <- Sys.time()
end_time - start_time

# Plot smoothed DMR methylation -------------------------------------------
 
glue::glue("Annotating DMRs and plotting smoothed values...")
pData <- pData(bs.filtered.bsseq)
if(length(levels(pData[,testCovariate])) == 2){
  pData$col <- NULL
  pData$col[pData[,testCovariate] == levels(pData[,testCovariate])[1]] <- "#3366CC"
  pData$col[pData[,testCovariate] == levels(pData[,testCovariate])[2]] <- "#FF3366"
  pData(bs.filtered.bsseq) <- pData
}
 
pdf("DMRs.pdf", height = 4, width = 8)
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

bs.filtered.bsseq %>%
  globalStats() %>%
  write.xlsx("smoothed_globalStats.xlsx") 

# PCAs of 20kb windows and CpG islands ------------------------------------

glue::glue("Creating and plotting PCA of 20 kb windows from {genome}")
goi %>%
  GenomeInfoDb::seqlengths() %>%
  GenomicRanges::tileGenome(tilewidth = 2e4,
                            cut.last.tile.in.chrom = TRUE) %>%
  GenomeInfoDb::keepStandardChromosomes(pruning.mode = "coarse") %>%
  cbind(., data.frame(
    bsseq::getMeth(BSseq = bs.filtered.bsseq,
                                     regions = .,
                                     type = "smooth",
                                     what = "perRegion"),
                      check.names = FALSE)
        ) %>%
  dplyr::select(-seqnames, -start, -end, -width, -strand) %>% 
  na.omit() %>%
  as.matrix() %>%
  t() %>% 
  DMRichR::PCA(group = bs.filtered.bsseq %>%
                 pData() %>%
                 as.tibble() %>%
                 pull(!!testCovariate),
               title = "Smoothed 20 Kb CpG Windows with CpG Islands") %>%
  ggsave("Smoothed 20 Kb CpG Windows with CpG Islands.pdf",
         plot = .,
         device = NULL)
 
if(genome == "hg38" | genome == "mm10" | genome == "rn6"){
  glue::glue("Creating and plotting PCA of CpG islands from {genome}")
  annotatr::build_annotations(genome = genome,
                           annotations = paste(genome,"_cpg_islands", sep = "")) %>% 
    GenomeInfoDb::keepStandardChromosomes(pruning.mode = "coarse") %>% 
    cbind(., data.frame(
      bsseq::getMeth(BSseq = bs.filtered.bsseq,
                                       regions = .,
                                       type = "smooth",
                                       what = "perRegion"),
                        check.names = FALSE)
          ) %>% 
    dplyr::select(-seqnames, -start, -end, -width, -strand,
                  - id, -tx_id, -gene_id, -symbol, - type) %>% 
    na.omit() %>%
    as.matrix() %>%
    t() %>% 
    DMRichR::PCA(group = bs.filtered.bsseq %>%
                   pData() %>%
                   as.tibble() %>%
                   pull(!!testCovariate),
                 title = "Smoothed CpG Island Windows") %>%
    ggsave("Smoothed CpG Island Windows.pdf",
           plot = .,
           device = NULL)
}

# Heatmap -----------------------------------------------------------------

pdf("heatmap.pdf", height = 8.5, width = 11)
smoothHeatmap(regions = sigRegions,
              bsseq = bs.filtered.bsseq,
              groups = bs.filtered.bsseq %>% pData() %>% as.tibble() %>% pull(!!testCovariate),
              out = "sig_individual_smoothed_DMR_methylation.txt")
dev.off()

# Prepare files for enrichment analyses -----------------------------------

if(genome == "hg38" | genome == "mm10" | genome == "rn6"){
  cat("\n[DMRichR] Preparing files for annotations \t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
  
  glue::glue("Preparing DMRs for annotations...")
  external <- sigRegions %>%
    labelDirection()
  
  glue::glue("Preparing background regions for annotations...")
  external_bg <- regions %>%
    labelDirection()
  
  glue::glue("Preparing regions for external GAT analysis...")
  dir.create("GAT")
  
  external %>%
    GenomeInfoDb::as.data.frame() %>%
    dplyr::select(seqnames, start, end, direction) %>% 
    df2bed("GAT/DMRs.bed")
  
  external_bg %>%
    GenomeInfoDb::as.data.frame() %>%
    dplyr::select(seqnames, start, end) %>% 
    df2bed("GAT/background.bed")

  glue::glue("Preparing DMRs for external HOMER analysis...")
  dir.create("HOMER")
  
  external %>%
    GenomeInfoDb::as.data.frame() %>% 
    dplyr::filter(direction == "Hypermethylated") %>%
    dplyr::select(seqnames, start, end) %>%
    df2bed("HOMER/DMRs_hyper.bed")
  
  external %>%
    GenomeInfoDb::as.data.frame() %>% 
    dplyr::filter(direction == "Hypomethylated") %>%
    dplyr::select(seqnames, start, end) %>%
    df2bed("HOMER/DMRs_hypo.bed")
  
  external_bg %>%
    GenomeInfoDb::as.data.frame() %>%
    dplyr::select(seqnames, start, end) %>% 
    df2bed("HOMER/background.bed")

  # CpG annotations ---------------------------------------------------------

  cat("\n[DMRichR] Building CpG annotations \t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
  annotations <- build_annotations(genome = genome, annotations = paste(genome,"_cpgs", sep=""))
  annotations <- GenomeInfoDb::keepStandardChromosomes(annotations, pruning.mode = "coarse")

  glue::glue("Annotating DMRs...")
  dm_annotated_CpG <- annotate_regions(
    regions = external,
    annotations = annotations,
    ignore.strand = TRUE,
    quiet = FALSE)

  glue::glue("Annotating background regions...")
  background_annotated_CpG <- annotate_regions(
    regions = external_bg,
    annotations = annotations,
    ignore.strand = TRUE,
    quiet = FALSE)

  glue::glue("Saving files for GAT...")
  CpGs <- as.data.frame(annotations)
  CpGs <- CpGs[!grepl("_", CpGs$seqnames), ]
  table(CpGs$seqnames)
  df2bed(CpGs[, c(1:3,10)], paste("GAT/", genome, "CpG.bed", sep = ""))

  glue::glue("Preparing CpG annotation plot...")
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
  ggsave("CpG_annotations.pdf", plot = CpG_bar, device = NULL, width = 8.5, height = 11)

  # Gene annotations --------------------------------------------------------

  cat("\n[DMRichR] Building gene region annotations \t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
  annotations <- build_annotations(genome = genome, annotations = c(paste(genome,"_basicgenes", sep = ""),
                                                                    paste(genome,"_genes_intergenic", sep = ""),
                                                                    paste(genome,"_genes_intronexonboundaries", sep = ""),
                                                                    if(genome == "hg38" | genome == "mm10"){paste(genome,"_enhancers_fantom", sep = "")}))
  annotations <- GenomeInfoDb::keepStandardChromosomes(annotations, pruning.mode = "coarse")

  glue::glue("Saving files for GAT...")
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

  glue::glue("Annotating DMRs...")
  dm_annotated <- annotate_regions(
    regions = external,
    annotations = annotations,
    ignore.strand = TRUE,
    quiet = FALSE)

  glue::glue("Annotating background regions...")
  background_annotated <- annotate_regions(
    regions = external_bg,
    annotations = annotations,
    ignore.strand = TRUE,
    quiet = FALSE)

  glue::glue("Preparing CpG annotation plot...")
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
    scale_y_continuous(expand = c(0,0)) +
    theme_classic() +
    theme(axis.text = element_text(size = 25),
          axis.title = element_text(size = 25),
          strip.text = element_text(size = 25),
          legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave("generegion_annotations.pdf", plot = gene_bar, device = NULL, width = 8.5, height = 11)
}

# ChIPseeker --------------------------------------------------------------

cat("\n[DMRichR] Annotating DMRs with gene symbols \t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
peakAnno <- annotatePeak(sigRegions,
                         TxDb = TxDb,
                         annoDb = annoDb,
                         overlap = "all")

glue::glue("Annotating background regions with gene symbols...")
backgroundAnno <- annotatePeak(regions,
                               TxDb = TxDb,
                               annoDb = annoDb,
                               overlap = "all")

glue::glue("Upset Plot of genic features...")
pdf("Upset.pdf", onefile = FALSE, height = 8.5, width = 11)
upsetplot(peakAnno, vennpie = TRUE)
dev.off()

glue::glue("Saving gene annotations...")

# Excel
peakAnno %>%
  tidyDMRs() %>%
  write.xlsx(file = "DMRs_annotated.xlsx", sep= "")

backgroundAnno %>%
  tidyDMRs() %>%
  write.xlsx(file = "background_annotated.xlsx", sep= "")

# Html
peakAnno %>%
  DMReport()

# Manhattan and Q-Q plots -------------------------------------------------

backgroundAnno %>%
  manQQ()

# Gene ontology and pathway analyses  -------------------------------------

cat("\n[DMRichR] Performing gene ontology and pathway analyses \t", format(Sys.time(), "%d-%m-%Y %X"), "\n")

glue::glue("Running enrichR")
# Check available databases
#dbs <- listEnrichrDbs()
dbs <- c("GO_Biological_Process_2018",
         "GO_Cellular_Component_2018",
         "GO_Molecular_Function_2018",
         "KEGG_2016",
         "Panther_2016",
         "Reactome_2016",
         "RNA-Seq_Disease_Gene_and_Drug_Signatures_from_GEO")

enrichResults <- peakAnno %>%
  tidyDMRs() %>%
  dplyr::select(geneSymbol) %>%
  purrr::flatten() %>%
  enrichr(dbs)

enrichResults %>%
  GOplot(tool = "enrichR") %>%
  ggsave("enrichr_plot.pdf",
         plot = .,
         device = NULL,
         height = 8.5,
         width = 12)

enrichResults %>%
  write.xlsx(file = "enrichr.xlsx")


glue::glue("Running rGREAT")
if(genome == "hg38" | genome == "mm10"){
  
  GREATjob <- GREAT(sigRegions = sigRegions,
                    regions = regions,
                    genome = genome)
  
  glue::glue("Saving and plotting GREAT enrichment results...")
  
  GREATresults <- GREATjob %>%
    rGREAT::getEnrichmentTables(category = c("GO", "Pathway Data"))
  
  GREATresults %>% 
    GOplot(tool = "rGREAT") %>%
    ggsave("GREAT_plot.pdf",
           plot = .,
           device = NULL,
           height = 8.5,
           width = 12)
  
  GREATresults %>%
    write.xlsx(file = "GREAT_results.xlsx")
  
  glue::glue("Saving and plotting GREAT gene annotation data...")
  pdf("GREAT_gene_associations_graph.pdf",
      height = 8.5,
      width = 11)
  par(mfrow = c(1, 3))
  res <- plotRegionGeneAssociationGraphs(GREATjob)
  dev.off()
  write.csv(as.data.frame(res),
            file = "GREATannotations.csv",
            row.names = F)
}


glue::glue("Running GOfuncR")
GOfuncResults <- GOfuncR(sigRegions = sigRegions,
                         regions = regions,
                         genome = genome,
                         upstream = 5000,
                         downstream = 1000,
                         annoDb = annoDb,
                         TxDb = TxDb)

GOfuncResults %>% 
  GOplot(tool = "GOfuncR") %>% 
  ggsave("GOfuncR_plot.pdf",
       plot = .,
       device = NULL,
       height = 8.5,
       width = 12)

GOfuncResults %>%
  write.xlsx("GOfuncR.xlsx")

glue::glue("Saving RData...")
GO_env <- ls(all = TRUE)[!(ls(all = TRUE) %in% bismark_env) &
                           !(ls(all = TRUE) %in% DMRs_env) &
                           !(ls(all = TRUE) %in% bsseq_env)]
save(list = GO_env, file = "GO.RData")
#load("GO.RData")
                          
# Blocks ------------------------------------------------------------------

cat("\n[DMRichR] Testing for large blocks (PMDs/HMDs) \t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
start_time <- Sys.time()

blocks <- dmrseq(bs = bs.filtered,
                 cutoff = cutoff,
                 maxPerms = maxPerms,
                 testCovariate = testCovariate,
                 adjustCovariate = adjustCovariate,
                 matchCovariate = matchCovariate,
                 block = TRUE,
                 minInSpan = 500,
                 bpSpan = 5e4,
                 maxGapSmooth = 1e6,
                 maxGap = 5e3)

glue::glue("Selecting significant blocks...")

if(sum(blocks$qval < 0.05) == 0 & sum(blocks$pval < 0.05) != 0){
  sigBlocks <- blocks[blocks$pval < 0.05,]
}else if(sum(blocks$qval < 0.05) >= 1){
  sigBlocks <- blocks[blocks$qval < 0.05,]
}else if(sum(blocks$pval < 0.05) == 0 & length(blocks) != 0){
  print("No significant blocks detected")
}else if(length(blocks) == 0){
  glue::glue("No background blocks detected, workflow is complete")
  quit(save = "no", status = 0, runLast = FALSE)
}

glue::glue("{length(sigBlocks)} significant blocks of differential methylation in {length(blocks)} background blocks")

glue::glue("Exporting block and background information...")
gr2csv(blocks, "backgroundBlocks.csv")
gr2bed(blocks, "backgroundBlocks.bed")
if(sum(blocks$pval < 0.05) > 0){
  gr2csv(sigBlocks, "blocks.csv")
  gr2bed(sigBlocks, "blocks.bed")
}

if(sum(blocks$pval < 0.05) > 0){
  glue::glue("Annotating and plotting blocks...")
  pdf("Blocks.pdf", height = 7.50, width = 11.50)
  plotDMRs(bs.filtered,
           regions = sigBlocks,
           testCovariate = testCovariate,
           annoTrack = getAnnot(genome),
           qval = FALSE)
  dev.off()
}

glue::glue("Saving RData...")
blocks_env <- ls(all = TRUE)[!(ls(all = TRUE) %in% bismark_env) &
                               !(ls(all = TRUE) %in% DMRs_env) &
                               !(ls(all = TRUE) %in% bsseq_env) &
                               !(ls(all = TRUE) %in% GO_env)]
save(list = blocks_env, file = "Blocks.RData")
#load("Blocks.RData")

glue::glue("Blocks timing...")
end_time <- Sys.time()
end_time - start_time

# End ---------------------------------------------------------------------

cat("\n[DMRichR] Summary \t\t\t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")

glue::glue("{length(sigRegions)} Significant DMRs \\
           ({round(sum(sigRegions$stat > 0) / length(sigRegions), digits = 2)*100}% hypermethylated, \\
           {round(sum(sigRegions$stat < 0) / length(sigRegions), digits = 2)*100}% hypomethylated) \\
           in {length(regions)} background regions \\
           from {nrow(bs.filtered)} CpGs assayed at {coverage}x coverage")

glue::glue("{length(sigBlocks)} significant blocks of differential methylation in {length(blocks)} background blocks")

if(length(grep("genomecenter.ucdavis.edu", .libPaths())) == 0){sessionInfo()}
rm(list = ls())
glue::glue("Done...")
quit(save = "no", status = 0, runLast = FALSE)

