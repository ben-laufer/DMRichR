#!/usr/bin/env Rscript

# DM.R
# Ben Laufer

rm(list=ls())
options(scipen=999)

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
  stopifnot(suppressMessages(sapply(CRAN, require, character.only= TRUE)))
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
    BiocManager::install(new.packages, ask = FALSE, quiet = TRUE)
  }
  stopifnot(suppressMessages(sapply(packages, require, character.only = TRUE)))
}

#' cleanRanges
#' @description Remove junk contigs and the mitochondrial chromosome from a genomic ranges or bsseq object.
#' @param gr Genomic ranges or bsseq object
#' @return Genomic ranges or bsseq object with junk contigs and the mitochondrial chromosome removed
#' @export cleanRanges
cleanRanges <- function(gr = gr){
  cat("\n[DM.R] Removing junk contigs and mitochondrial DNA \t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
  stopifnot(is(gr, "BSseq") | is(gr, "GRanges"))
  print(length(seqlevels(gr)))
  gr <- keepStandardChromosomes(gr, pruning.mode = "coarse")
  gr <- dropSeqlevels(gr, "chrM", pruning.mode = "coarse")
  print(length(seqlevels(gr)))
  return(gr)
}

#' gr2csv
#' @description Save a genomic ranges object as a csv file
#' @param gr Genomic ranges or bsseq object
#' @param csv Name of the csv file in quotations
#' @return CSV file
#' @export gr2csv
gr2csv <- function(gr = gr, csv = csv){
  write.csv(as.data.frame(gr), file = csv, row.names = FALSE) 
}

#' gr2bed
#' @description Save a genomic ranges object as a basic bed file
#' @param gr Genomic ranges or bsseq object
#' @param bed Name of the bed file in quotations
#' @return Bed file
#' @export gr2bed
gr2bed <- function(gr = gr, bed = bed){
  write.table(as.data.frame(gr)[1:3], bed, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
}

#' getSmooth
#' @description Provides individual smoothed methylation values for genomic ranges objects using bsseq
#' @param bsseq Smoothed bsseq object
#' @param regions Genomic ranges object 
#' @param out Name of the text file in quotations
#' @return Genomic ranges object of individual smoothed methylation values and text file
#' @export getSmooth
getSmooth <- function(bsseq = bsseq,
                      regions = regions,
                      out = out){
  smoothed <- data.frame(getMeth(BSseq = bsseq, regions = regions, type = "smooth", what = "perRegion"))
  colnames(smoothed) <- names
  smoothed_table <- cbind(regions, smoothed)
  write.table(smoothed_table, out, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  return(smoothed_table)
}

#' smooth2txt
#' @description Save smoothed methylation values as a text file
#' @param df Data frame
#' @param txt Name of the text file in quotations
#' @return Text file
#' @export smooth2txt
smooth2txt <- function(df = df, txt = txt){
  write.table(df, txt, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}

#' PCA
#' @description Provides individual smoothed methylation values for genomic ranges objects using bsseq
#' @param matrix Matrix of transposed individual methylation values
#' @param title Character string of title for plot and pdf
#' @return PCA plot
#' @export PCA
PCA <- function(matrix = matrix,
                title = title){
  
  cat("\n[DM.R] Performing PCA \t\t\t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
  data.pca <- prcomp(matrix, center = TRUE, scale. = TRUE) 
  plot(data.pca, type = "l")
  print(summary(data.pca))
  
  cat("\n[DM.R] Plotting PCA \t\t\t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
  PCA <- ggbiplot(data.pca,
                  obs.scale = 1,
                  var.scale = 1,
                  groups = group,
                  ellipse = TRUE,
                  circle = FALSE,
                  var.axes = FALSE,
                  choices = 1:2) +
    scale_color_discrete(name = '') +
    theme_bw(base_size = 25) +
    geom_point(aes(colour = group), size = 4) +
    theme(legend.direction = 'vertical',
          legend.position = c(0.125, 0.1), # Change legend position
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 18),
          panel.grid.major = element_blank(), 
          panel.border = element_rect(color = "black", size = 1.25),
          axis.ticks = element_line(size = 1.25), 
          legend.key = element_blank(),
          panel.grid.minor = element_blank()) +
    guides(col=guide_legend(ncol=2)) +
    ggtitle(title) + # Change title
    theme(plot.title = element_text(hjust = 0.5)) 
  ggsave(paste(title,".pdf", sep = ""), plot = PCA, device = NULL)
  return(PCA)
}

#' df2bed
#' @description Save a dataframe as a basic bed file
#' @param df Data frame
#' @param bed Name of the bed file in quotations
#' @return Bed file
#' @export df2bed
df2bed <-function(df = df, bed = bed){
  write.table(df, bed, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
}

#https://stackoverflow.com/questions/8197559/emulate-ggplot2-default-color-palette
#' gg_color_hue
#' @description Generate ggplot2 style colors
#' @param n Number of samples
#' @return Character string of colors
#' @export gg_color_hue
gg_color_hue <- function(n = n){
  cat("\n[DM.R] Preparing Colors \t\t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

# Install and update ------------------------------------------------------

cat("\n[DM.R] Installing and updating pacakges \t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
packageManage()
packageLoad(c("tidyverse", "dmrseq", "annotatr", "rGREAT", "enrichR", "ChIPseeker", "BiocParallel", "ggbiplot",
              "liftOver", "openxlsx", "CMplot", "optparse", "devtools", "gplots", "RColorBrewer", "nlme", "lsmeans"))
suppressWarnings(BiocManager::valid(fix = TRUE, update = TRUE, ask = FALSE))

# Global variables --------------------------------------------------------

cat("\n[DM.R] Processing arguments from script \t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")

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
  make_option(c("-c", "--cores"), type = "integer", default = 2,
              help = "Choose number of cores [default = %default]")
)
opt <- parse_args(OptionParser(option_list=option_list))

cat("\n[DM.R] Assigning arguments to global variables \t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
# Check for requirements
stopifnot(!is.null(opt$genome))
stopifnot(!is.null(opt$testCovariate))
# Assign
genome <- as.character(opt$genome)
coverage <- as.numeric(opt$coverage)
testCovariate <- as.character(opt$testCovariate)
if(!is.null(opt$adjustCovariate)){
  adjustCovariate <- as.character(opt$adjustCovariate)
}else if(is.null(opt$adjustCovariate)){
  adjustCovariate <- opt$adjustCovariate
}
if(!is.null(opt$matchCovariate)){
  matchCovariate <- as.character(opt$matchCovariate)
}else if(is.null(opt$matchCovariate)){
  matchCovariate <- opt$matchCovariate
}
cores <- as.numeric(opt$cores)

# Setup Annotation Databases ----------------------------------------------

cat("\n[DM.R] Selecting annotation databases \t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")

if(genome == "hg38"){
  goi <- BSgenome.Hsapiens.UCSC.hg38; TxDb <- TxDb.Hsapiens.UCSC.hg38.knownGene; annoDb <- "org.Hs.eg.db"
  packages <- c("TxDb.Hsapiens.UCSC.hg38.knownGene", "org.Hs.eg.db", "BSgenome.Hsapiens.UCSC.hg38")
}else if(genome == "mm10"){
  goi <- BSgenome.Mmusculus.UCSC.mm10; TxDb <- TxDb.Mmusculus.UCSC.mm10.knownGene; annoDb <- "org.Mm.eg.db"
  packages <- c("TxDb.Mmusculus.UCSC.mm10.knownGene", "org.Mm.eg.db", "BSgenome.Mmusculus.UCSC.mm10")
}else if(genome == "rheMac8"){
  goi <- BSgenome.Mmulatta.UCSC.rheMac8; TxDb <- TxDb.Mmulatta.UCSC.rheMac8.refGene; annoDb <- "org.Mmu.eg.db"
  packages <- c("TxDb.Mmulatta.UCSC.rheMac8.refGene", "org.Mmu.eg.db", "BSgenome.Mmulatta.UCSC.rheMac8")
}else if(genome == "rn6"){
  goi <- BSgenome.Rnorvegicus.UCSC.rn6; TxDb <- TxDb.Rnorvegicus.UCSC.rn6.refGene; annoDb <- "org.Rn.eg.db"
  packages <- c("TxDb.Rnorvegicus.UCSC.rn6.refGene", "org.Rn.eg.db", "BSgenome.Rnorvegicus.UCSC.rn6")
}else{
  stop(paste(genome, "is not suppourted, please choose either hg38, mm10, rheMac8, or rn6 [Case Sensitive]"))
}
packageLoad(packages)

# Load and process samples ------------------------------------------------

cat("\n[DM.R] Loading Bismark cytosine reports \t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
cov <- list.files(path=getwd(), pattern="*.txt.gz") 
names <- gsub( "_.*$","", cov)

bs <- read.bismark(files = cov,
                   sampleNames = names,
                   rmZeroCov = TRUE,
                   strandCollapse = TRUE,
                   fileType = "cytosineReport",
                   verbose = TRUE,
                   mc.cores = cores)

cat("\n[DM.R] Assigning sample metadata \t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
meta <- read.csv("sample_info.csv", header = TRUE)
meta <- meta[order(match(meta[,1],names)),]
stopifnot(sampleNames(bs) == meta$Name)
pData(bs) <- cbind(pData(bs), meta[2:length(meta)])
pData(bs)

cat("\n[DM.R] Filtering CpGs for coverage \t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
bs <- cleanRanges(bs)
bs
head(getCoverage(bs, type = "Cov"))
sample.idx <- which(pData(bs)[[testCovariate]] %in% levels(pData(bs)[[testCovariate]]))
loci.idx <- which(DelayedMatrixStats::rowSums2(getCoverage(bs, type="Cov") >= coverage) >= length(sample.idx))
bs.filtered <- bs[loci.idx, sample.idx]
bs.filtered
head(getCoverage(bs.filtered, type = "Cov"))

cat("\n[DM.R] Saving Rdata \t\t\t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
bismark_env <- ls(all = TRUE)
save(list = bismark_env, file = "bismark.RData") 
#load("bismark.RData")

# Distribtuion plots ------------------------------------------------------

cat("\n[DM.R] Plotting Empirical Distribution of CpGs \t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
pdf("Filtered_CpG_Methylation_Distributions.pdf", height = 7.50, width = 11.50)
plotEmpiricalDistribution(bs.filtered, testCovariate = testCovariate)
plotEmpiricalDistribution(bs.filtered, testCovariate = testCovariate, type = "Cov", bySample = TRUE)
dev.off()

# DMRs --------------------------------------------------------------------

cat("\n[DM.R] Testing for DMRs with dmrseq \t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
# Reproducible permutations (change and record seed for different datasets to avoid any potential random bias)
set.seed(1)
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

cat("\n[DM.R] Selecting significant DMRs \t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
if(sum(regions$qval < 0.05) < 100){
  sigRegions <- regions[regions$pval < 0.05,]
}else if(sum(regions$qval < 0.05) >= 100){
  sigRegions <- regions[regions$qval < 0.05,]  
}else if(sum(regions$pval < 0.05) == 0){
  stop("No significant DMRs detected")  
  }
cat(paste(round(sum(sigRegions$stat > 0) / length(sigRegions), digits = 2)*100, "% of DMRs are hypermethylated", sep =""))

cat("\n[DM.R] Plotting DMR pie chart \t\t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
pie <- (table(sigRegions$stat < 0))
names(pie) <- c("Hypermethylated", "Hypomethylated")
pdf("HypervsHypo_pie.pdf", height = 8.5, width = 11)
pie(pie,
    labels = c(paste(pie[1], "Hypermethylated", sep = " "), paste(pie[2], "Hypomethylated", sep = " ")),
    col = c("Red", "Blue"))
dev.off()

cat("\n[DM.R] Extracing raw differnces for DMRs \t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
rawDiff <- meanDiff(bs.filtered,
                    dmrs = regions,
                    testCovariate = testCovariate)
sigRawDiff <- meanDiff(bs.filtered,
                       dmrs = sigRegions,
                       testCovariate = testCovariate)
regions$RawDiff <- rawDiff
sigRegions$RawDiff <- sigRawDiff

cat("\n[DM.R] Exporting DMR and background region information \t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
gr2csv(regions, "backgroundRegions.csv")
gr2csv(sigRegions, "DMRs.csv")
gr2bed(regions, "backgroundRegions.bed")
gr2bed(sigRegions, "DMRs.bed")

cat("\n[DM.R] Annotating and plotting DMRs \t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
pdf("DMRs.pdf", height = 7.50, width = 11.50)
annoTrack <- getAnnot(genome)
plotDMRs(bs.filtered, regions = sigRegions, testCovariate = testCovariate, annoTrack = annoTrack, qval = F)
dev.off()

cat("\n[DM.R] Saving Rdata \t\t\t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
DMRs_env <- ls(all = TRUE)[!(ls(all = TRUE) %in% bismark_env)]
save(list = DMRs_env, file = "DMRs.RData")
#load("DMRs.RData")

# Individual smoothed values ----------------------------------------------

cat("\n[DM.R] Smoothing individual methylation values \t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
bs.filtered.bsseq <- BSmooth(bs.filtered, mc.cores = cores, verbose = TRUE)
bs.filtered.bsseq

cat("\n[DM.R] Extracting values for WGCNA \t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
indiv_smoothed_table <- getSmooth(bsseq = bs.filtered.bsseq,
                                  regions = regions,
                                  out = "background_region_individual_smoothed_methylation.txt")

cat("\n[DM.R] Saving Rdata \t\t\t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
bsseq_env <- ls(all = TRUE)[!(ls(all = TRUE) %in% bismark_env) &
                              !(ls(all = TRUE) %in% DMRs_env)]
save(list = bsseq_env, file = "bsseq.RData") 
#load("bsseq.RData")

# Global methylation ------------------------------------------------------
# Warning! not finished
if(length(adjustCovariate) == 1){
  cat("\n[DM.R] Extracting global methylation \t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
  # Global for all chromosomes
  global <- data.frame(DelayedMatrixStats::colMeans2(getMeth(BSseq = bs.filtered.bsseq, type = "smooth", what = "perBase")))
  global$sample <- names
  names(global) <- c("CpG_Avg", "sample")
  global <- as.tibble(cbind(global, data.frame(pData(bs.filtered.bsseq))), rownames = NULL) %>% 
    select(sample,
           CpG_Avg,
           testCovariate,
           adjustCovariate,
           matchCovariate) %>%
    rename(sample ="sample",
           testCovariate = !!testCovariate,
           adjustCovariate = !!adjustCovariate,
           matchCovariate = !!matchCovariate)
  global       
  
  cat("\n[DM.R] Fitting linear model \t\t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n") 
  if(length(levels(global$matchCovariate)) == 1){
    lmefit <- lme(CpG_Avg ~ testCovariate + adjustCovariate, ~1|sample, data = global)
  }else if(length(levels(global$matchCovariate)) > 1){
    lmefit <-lme(CpG_Avg ~ testCovariate + matchCovariate + adjustCovariate, ~1|sample, data = global)
  }
  
  cat("\n[DM.R] ANOVA and post-hoc comparisons \t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n") 
  anova <- lmefit %>%
    anova(type = "marginal") %>% 
    rownames_to_column()  %>%
    as.tibble()  
  postHoc <- lmefit %>%
    ref.grid() %>%
    lsmeans(~testCovariate) %>%
    pairs() %>% 
    summary() %>%
    as.tibble()
  write.xlsx(list("anova" = anova,
                  "postHoc" = postHoc,
                  "lme" = broom::augment(lmefit)),
             "smoothed_global_methylation_stats.xlsx")
}

# Chromosomal methylation -------------------------------------------------
# Warning! not finished
if(length(adjustCovariate) == 1){
  cat("\n[DM.R] Extracting chromosomal methylation \t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
  grl <- split(bs.filtered.bsseq, seqnames(bs.filtered.bsseq))
  global_chr <- matrix(ncol = length((seqlevels(grl))), nrow = 1)
  for(i in seq_along(seqlevels(grl))){
    global_chr[i] <- data.frame(DelayedMatrixStats::colMeans2(getMeth(BSseq = grl[[i]], type = "smooth", what = "perBase")))
    names(global_chr)[i] <- seqlevels(grl)[i]
  }
  global_chr$sample <- names
  global_chr <- as.tibble(cbind(global_chr, data.frame(pData(bs.filtered.bsseq))), rownames = NULL) %>%
    rename(sample ="sample",
           testCovariate = !!testCovariate,
           adjustCovariate = !!adjustCovariate,
           matchCovariate = !!matchCovariate) %>%
    select(sample,
           testCovariate,
           adjustCovariate,
           matchCovariate,
           contains("chr")) %>%
    gather(key = chromosome,
           value = CpG_Avg,
           -sample,
           -testCovariate,
           -adjustCovariate,
           -matchCovariate)
  global_chr
  
  cat("\n[DM.R] Fitting linear model \t\t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n") 
  if(length(levels(global_chr$matchCovariate)) == 1){
    lmefit <-lme(CpG_Avg ~ testCovariate + adjustCovariate + chromosome + testCovariate*chromosome + testCovariate*adjustCovariate, ~1|sample, data = global_chr)
  }else if(length(levels(global_chr$matchCovariate)) > 1){
    lmefit <-lme(CpG_Avg ~ testCovariate + adjustCovariate + matchCovariate + chromosome + testCovariate*chromosome + testCovariate*adjustCovariate + testCovariate*matchCovariate, ~1|sample, data = global_chr)
  }
  
  cat("\n[DM.R] ANOVA and post-hoc comparisons \t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
  anova <- lmefit %>%
    anova(type = "marginal") %>% 
    rownames_to_column() %>%
    as.tibble()  
  postHoc <- lmefit %>%
    ref.grid() %>%
    lsmeans(~testCovariate|chromosome) %>%
    pairs() %>%
    summary() %>%
    as.tibble() %>%
    mutate(fdr = p.adjust(p.value, method = 'fdr'))
  write.xlsx(list("anova" = anova,
                  "postHoc" = postHoc,
                  "lme" = broom::augment(lmefit)),
             "smoothed_global_chromosomal_methylation_stats.xlsx")
}

# PCA of 20 kb windows with CGi -------------------------------------------

cat("\n[DM.R] Creating 20 kb windows \t\t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
chrSizes <- seqlengths(goi)
windows <- tileGenome(chrSizes,
                      tilewidth = 2e4,
                      cut.last.tile.in.chrom = TRUE)
windows <- cleanRanges(windows)
windows

cat("\n[DM.R] Extracting values for 20 kb windows \t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
windows_smoothed_table <- getSmooth(bsseq = bs.filtered.bsseq,
                                    regions = windows,
                                    out = "20kb_smoothed_windows.txt")

cat("\n[DM.R] Tidying 20 kb window data \t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n") 
meth_reorder <- na.omit(windows_smoothed_table[,c(6:length(windows_smoothed_table))])
data <- t(as.matrix(meth_reorder))
# Fix for columns of no variance
#data2 <- as.data.frame(data)
#data3 <- data2[,apply(data2, 2, var, na.rm=TRUE) != 0]
# Titles
stopifnot(sampleNames(bs.filtered.bsseq) == colnames(meth_reorder))
group <- as.tibble(pData(bs.filtered.bsseq)) %>% pull(!!testCovariate)

cat("\n[DM.R] 20 kb window PCA \t\t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
PCA(data, "Smoothed 20 Kb CpG Windows with CpG Islands")

# PCA of CGi windows ------------------------------------------------------

if(genome == "hg38" | genome == "mm10" | genome == "rn6"){
  cat("\n[DM.R] Creating CGi windows \t\t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
  annots <- paste(genome,"_cpg_islands", sep = "")
  CGi <- build_annotations(genome = genome, annotations = annots)
  CGi <- cleanRanges(CGi)
  CGi
  
  cat("\n[DM.R] Extracting values for CGi windows \t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
  CGi_smoothed_table <- getSmooth(bsseq = bs.filtered.bsseq,
                                  regions = CGi,
                                  out = "CGi_smoothed_windows.txt")
  
  cat("\n[DM.R] Tidying CGi window data \t\t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
  meth_reorder <- na.omit(CGi_smoothed_table[,c(11:length(CGi_smoothed_table))])
  data <- t(as.matrix(meth_reorder))
  # Fix for columns of no variance
  #data2 <- as.data.frame(data)
  #data3 <- data2[,apply(data2, 2, var, na.rm=TRUE) != 0]
  # Titles
  stopifnot(sampleNames(bs.filtered.bsseq) == colnames(meth_reorder))
  group <- as.tibble(pData(bs.filtered.bsseq)) %>% pull(!!testCovariate)
  
  cat("\n[DM.R] CGi window PCA \t\t\t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n") 
  PCA(data, "Smoothed CpG Island Windows")
}

# Heatmap -----------------------------------------------------------------

#' smoothHeatmap
#' @description Plot a heatmap of invdividual smoothed methylation value z scores for significant DMRs
#' @param sigRegions class bsseq object of significant DMRs
#' @return Saves a pdf image of the heatmap
#' @export smoothHeatmap
smoothHeatmap <- function(sigRegions = sigRegions){
  cat("\n[DM.R] Extracting values for DMR heatmap \t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
  sig_indiv_smoothed_table <- getSmooth(bsseq = bs.filtered.bsseq,
                                        regions = sigRegions,
                                        out = "sig_individual_smoothed_DMR_methylation.txt")
  
  cat("\n[DM.R] Tidying for heatmap of HCA \t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
  # Load smoothed values
  matrix <- as.matrix(sig_indiv_smoothed)
  # Convert to Percent
  matrix <- matrix[,]*100
  # Subtract the mean methylation for each row/DMR
  data <- sweep(matrix, 1, rowMeans(matrix)) 
  # Tidy
  data <- as.matrix(data)
  colnames(data) <- as.tibble(pData(bs.filtered.bsseq)) %>% pull(!!testCovariate)
  
  cat("\n[DM.R] Plotting heatmap of HCA \t\t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
  pdf("heatmap.pdf", height = 8.5, width = 11)
  heatmap.2(data,
            Rowv= as.dendrogram(hclust(dist(data))),
            scale = c("row"),
            Colv = TRUE,
            col = rev(brewer.pal(11, name = "RdBu")),
            margins = c(10,10),
            trace = "none",
            main = paste(nrow(sig_indiv_smoothed),"Differentially Methylated Regions", sep = " "),
            labRow = NA,
            srtCol = 60,
            keysize = 0.85,
            key.par = list(cex=0.5),
            key.xlab = "Z-score(% mCG/CG - mean)",
            key.ylab = "Frequency",
            key.title = ""
  )
  dev.off()
}
smoothHeatmap(sigRegions)

# Prepare files for enrichment analyses -----------------------------------

if(genome == "hg38" | genome == "mm10" | genome == "rn6"){
  cat("\n[DM.R] Preparing DMRs files for annotations \t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
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
  write.table(externalOut[,c(1:3,16)], "GAT/DMRs.bed", quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
  
  cat("\n[DM.R] Preparing DMRs for HOMER \t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
  dir.create("HOMER")
  write.table((external[,c(1:3)])[external$direction == "Hypermethylated",], "HOMER/DMRs_hyper.bed", quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
  write.table((external[,c(1:3)])[external$direction == "Hypomethylated",], "HOMER/DMRs_hypo.bed", quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
  
  cat("\n[DM.R] Preparing background regions for annotations \t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
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
  
  cat("\n[DM.R] Building CpG annotations \t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
  annots <- paste(genome,"_cpgs", sep="")
  annotations <- build_annotations(genome = genome, annotations = annots)
  annotations <- cleanRanges(annotations)
  
  cat("\n[DM.R] Annotating DMRs \t\t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
  dm_annotated_CpG <- annotate_regions(
    regions = external,
    annotations = annotations,
    ignore.strand = TRUE,
    quiet = FALSE)
  
  cat("\n[DM.R] Annotating background regions \t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
  background_annotated_CpG <- annotate_regions(
    regions = external_bg,
    annotations = annotations,
    ignore.strand = TRUE,
    quiet = FALSE)
  
  cat("\n[DM.R] Saving files for GAT \t\t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
  CpGs <- as.data.frame(annotations)
  CpGs <- CpGs[!grepl("_", CpGs$seqnames), ]
  table(CpGs$seqnames)
  df2bed(CpGs[, c(1:3,10)], paste("GAT/", genome, "CpG.bed", sep = ""))
  
  cat("\n[DM.R] Preparing CpG annotation plot \t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
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
  
  cat("\n[DM.R] Building gene region annotations \t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
  annots <- c(paste(genome,"_basicgenes", sep = ""),
             paste(genome,"_genes_intergenic", sep = ""),
             paste(genome,"_genes_intronexonboundaries", sep = ""),
             if(genome == "hg38" | genome == "mm10"){paste(genome,"_enhancers_fantom", sep = "")})
  annotations <- build_annotations(genome = genome, annotations = annots)
  annotations <- cleanRanges(annotations)
  
  cat("\n[DM.R] Saving files for GAT \t\t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
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
  
  cat("\n[DM.R] Annotating DMRs \t\t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
  dm_annotated <- annotate_regions(
    regions = external,
    annotations = annotations,
    ignore.strand = TRUE,
    quiet = FALSE)
  
  cat("\n[DM.R] Annotating background regions \t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
  background_annotated <- annotate_regions(
    regions = external_bg,
    annotations = annotations,
    ignore.strand = TRUE,
    quiet = FALSE)
  
  cat("\n[DM.R] Preparing CpG annotation plot \t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
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
  cat("\n[DM.R] Obtaining liftOver information for GREAT \t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
  path <- system.file(package="liftOver", "extdata", "hg38ToHg19.over.chain")
  ch <- import.chain(path)
  ch
  
  cat("\n[DM.R] liftOver DMRs \t\t\t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
  seqlevelsStyle(sigRegions) <- "UCSC" 
  sigRegions_liftOver <- liftOver(sigRegions, ch)
  class(sigRegions_liftOver)
  sigRegions_liftOver <- unlist(sigRegions_liftOver)
  length(sigRegions) - length(sigRegions_liftOver)
  
  cat("\n[DM.R] liftOver background regions \t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
  seqlevelsStyle(regions) <- "UCSC" 
  regions_liftOver <- liftOver(regions, ch)
  class(regions_liftOver)
  regions_liftOver <- unlist(regions_liftOver)
  length(regions) - length(regions_liftOver)
}

if(genome == "hg38" | genome == "mm10"){
  cat("\n[DM.R] Submitting to GREAT\ \t\t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
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
  
  cat("\n[DM.R] Saving GREAT enrichment results \t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
  write.xlsx(tb, file = "GREAT_results.xlsx", sep="")
  
  cat("\n[DM.R] Plotting GREAT results \t\t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
  pdf("GREAT.pdf", height = 8.5, width = 11)
  par(mfrow = c(1, 3))
  res <- plotRegionGeneAssociationGraphs(job)
  dev.off()
  
  cat("\n[DM.R] Saving GREAT annotations \t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
  res
  write.csv(as.data.frame(res), file="GREATannotations.csv", row.names = F)
}

# Chipseeker --------------------------------------------------------------

cat("\n[DM.R] Annotating DMRs with gene symbols \t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
peakAnno <- annotatePeak(sigRegions,
                         TxDb = TxDb,
                         annoDb = annoDb,
                         overlap = "all")

backgroundAnno <- annotatePeak(regions,
                               TxDb = TxDb,
                               annoDb = annoDb,
                               overlap = "all")

cat("\n[DM.R] Upset Plot \t\t\t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
pdf("Upset.pdf", onefile = FALSE, height = 8.5, width = 11)
upsetplot(peakAnno, vennpie = TRUE)
dev.off()

cat("\n[DM.R] Saving gene annotations \t\t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
annotations <- as.data.frame(peakAnno)
annotations[2] <- annotations[2] - 1
annotations <- annotations[-c(5,17:19)]
write.xlsx(annotations, file = "DMRs_annotated.xlsx", sep= "")

background_annotations <- as.data.frame(backgroundAnno)
background_annotations[2] <- background_annotations[2] - 1
background_annotations <- background_annotations[-c(5,17:19)]
write.xlsx(background_annotations, file = "background_annotated.xlsx", sep= "")

# CMplot ------------------------------------------------------------------

cat("\n[DM.R] Tidying for Manhattan and QQ plots \t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
Manhattan <- as.data.frame(sort(as.GRanges(backgroundAnno), ignore.strand=TRUE))[c("SYMBOL","seqnames", "start", "pval")]
Manhattan$seqnames <- substring(Manhattan$seqnames, 4)
cols = gg_color_hue(2)

cat("\n[DM.R] Generating Manhattan and QQ plots \t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
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

cat("\n[DM.R] Running Enrichr \t\t\t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
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

cat("\n[DM.R] Saving RData \t\t\t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
GO_env <- ls(all = TRUE)[!(ls(all = TRUE) %in% bismark_env) &
                           !(ls(all = TRUE) %in% DMRs_env) &
                           !(ls(all = TRUE) %in% bsseq_env)]
save(list = GO_env, file = "GO.RData") 
#load("GO.RData")

# Blocks ------------------------------------------------------------------

cat("\n[DM.R] Testing for large blocks (PMDs/HMDs) \t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")

register(MulticoreParam(1))
blocks <- dmrseq(bs=bs.filtered,
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

cat("\n[DM.R] Selecting significant blocks \t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")

if(sum(blocks$qval < 0.05) == 0){
  sigBlocks <- blocks[blocks$pval < 0.05,]
}else if(sum(blocks$qval < 0.05) >= 1){
  sigBlocks <- blocks[blocks$qval < 0.05,]  
}else if(sum(blocks$pval < 0.05) == 0){
  warning("No significant DMRs detected")  
}

cat("\n[DM.R] Exporting block and background information \t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
gr2csv(blocks, "backgroundBlocks.csv")
gr2csv(sigBlocks, "blocks.csv")
gr2bed(blocks, "backgroundBlocks.bed")
gr2bed(sigBlocks, "blocks.bed")

cat("\n[DM.R] Annotating and plotting blocks \t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
pdf("Blocks.pdf", height = 7.50, width = 11.50)
annoTrack <- getAnnot(genome)
plotDMRs(bs.filtered, regions = sigBlocks, testCovariate = testCovariate, annoTrack = annoTrack, qval = FALSE)
dev.off()

cat("\n[DM.R] Saving RData \t\t\t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
blocks_env <- ls(all = TRUE)[!(ls(all = TRUE) %in% bismark_env) &
                               !(ls(all = TRUE) %in% DMRs_env) &
                               !(ls(all = TRUE) %in% bsseq_env) &
                               !(ls(all = TRUE) %in% GO_env)]
save(list = blocks_env, file = "Blocks.RData") 
#load("Blocks.RData")

# End ---------------------------------------------------------------------

cat("\n[DM.R] Finishing \t\t\t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
sessionInfo()
rm(list = ls())
cat("\n[DM.R] Done \t\t\t\t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
