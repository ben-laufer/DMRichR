#!/usr/bin/env Rscript

# DM.R
# Ben Laufer

# Version: 0.99.2
# Last update: September 22nd 2018

# Install and Update ------------------------------------------------------
rm(list=ls())
options(scipen=999)

cat("\n[DM.R] Installing and updating pacakges \t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")

CRAN <- c("BiocManager", "remotes")
new.packages <- CRAN[!(CRAN %in% installed.packages()[,"Package"])]
if(length(new.packages)>0){
  install.packages(new.packages, repos ="https://cloud.r-project.org")
  stopifnot(suppressMessages(sapply(CRAN, require, character.only= TRUE)))}

packages <- c("gplots","RColorBrewer","openxlsx","CMplot","stringr","ggplot2","plyr","devtools","dmrseq",
              "BiocParallel","annotatr","liftOver","rGREAT","enrichR","ChIPseeker",
              "TxDb.Hsapiens.UCSC.hg38.knownGene","org.Hs.eg.db","TxDb.Mmulatta.UCSC.rheMac8.refGene","org.Mmu.eg.db",
              "TxDb.Mmusculus.UCSC.mm10.knownGene","org.Mm.eg.db","TxDb.Rnorvegicus.UCSC.rn6.refGene","org.Rn.eg.db",
              "BSgenome.Hsapiens.UCSC.hg38","BSgenome.Mmulatta.UCSC.rheMac8","BSgenome.Mmusculus.UCSC.mm10",
              "BSgenome.Rnorvegicus.UCSC.rn6","vqv/ggbiplot", "optparse")

new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new.packages)>0){
  library("BiocManager")
  BiocManager::install(new.packages, ask = FALSE)
  packages <- gsub("vqv/", "", packages)
  stopifnot(suppressMessages(sapply(packages, require, character.only= TRUE)))}

suppressWarnings(BiocManager::valid(fix=TRUE, ask=FALSE))

# Global variables --------------------------------------------------------

cat("\n[DM.R] Processing arguments from script \t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")

option_list <- list( 
  make_option(c("-g", "--genome"), type = "character", default = NULL,
              help = "Choose a genome (hg38, mm10, rn6, rheMac8) [required]"),
  make_option(c("-x", "--coverage"), type = "integer", default = 1,
              help = "Choose a CpG coverage cutoff [default = 1]"),
  make_option(c("-t", "--testCovariate"), type = "character", default = NULL,
              help = "Choose a test covariate [required]"),
  make_option(c("-a", "--adjustCovariate"), type = "character", default = NULL,
              help = "Choose covariates to directly adjust [default = NULL]"),
  make_option(c("-m", "--matchCovariate"), type = "character", default = NULL,
              help = "Choose covariate to balance permutations [default = NULL]"),
  make_option(c("-c", "--cores"), type = "integer", default = 2,
              help = "Choose number of cores [default = 2]")
)

opt <- parse_args(OptionParser(option_list=option_list))

stopifnot(!is.null(opt$genome))
stopifnot(!is.null(opt$testCovariate))

genome <- as.character(opt$genome)
coverage <- as.numeric(opt$coverage)
testCovariate <- as.character(opt$testCovariate)
adjustCovariate <- as.character(opt$adjustCovariate)
matchCovariate <- as.character(opt$matchCovariate)
cores <- as.numeric(opt$cores)

# Setup Annotation Databases ----------------------------------------------

cat("\n[DM.R] Selecting annotation databases \t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")

if(genome == "hg38"){
  goi <- BSgenome.Hsapiens.UCSC.hg38; TxDb <- TxDb.Hsapiens.UCSC.hg38.knownGene; annoDb <- "org.Hs.eg.db"
}else if(genome == "mm10"){
  goi <- BSgenome.Mmusculus.UCSC.mm10; TxDb <- TxDb.Mmusculus.UCSC.mm10.knownGene; annoDb <- "org.Mm.eg.db"
}else if(genome == "rheMac8"){
  goi <- BSgenome.Mmulatta.UCSC.rheMac8; TxDb <- TxDb.Mmulatta.UCSC.rheMac8.refGene; annoDb <- "org.Mmu.eg.db"
}else if(genome == "rn6"){
  goi <- BSgenome.Rnorvegicus.UCSC.rn6; TxDb <- TxDb.Rnorvegicus.UCSC.rn6.refGene; annoDb <- "org.Rn.eg.db"
}else{
  stop(paste(genome, "is not suppourted, please choose either hg38, mm10, rheMac8, or rn6 [Case Sensitive"))
}

# Load and process samples ------------------------------------------------

cat("\n[DM.R] Loading Bismark cytosine reports \t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")

cov <- list.files(path=getwd(), pattern="*.txt.gz") 
names <- gsub( "_.*$","", cov)

bs<- read.bismark(files = cov,
                  sampleNames = names,
                  rmZeroCov = FALSE,
                  strandCollapse = TRUE,
                  fileType = "cytosineReport",
                  verbose = TRUE,
                  mc.cores = cores)

cat("\n[DM.R] Assigning sample metadata \t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
meta <- read.csv("sample_info.csv", header = TRUE)
meta <- meta[order(match(meta[,1],names)),]
stopifnot(sampleNames(bs) == meta$Name)
#tryCatch(stopifnot(sampleNames(bs) == meta$Name, error=stop("Samples names do not match between files and design matrix")))
pData(bs) <- cbind(pData(bs), meta[2:length(meta)])
pData(bs)

cat("\n[DM.R] Removing junk contigs and mitochondrial DNA \t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
length(seqlevels(bs))
bs <- keepStandardChromosomes(bs, pruning.mode = "coarse")
bs <- dropSeqlevels(bs, "chrM", pruning.mode = "coarse")
length(seqlevels(bs))

cat("\n[DM.R] Filtering CpGs for coverage \t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
bs
head(getCoverage(bs, type="Cov"))
sample.idx <- which(pData(bs)[[testCovariate]] %in% levels(pData(bs)[[testCovariate]]))
loci.idx <- which(DelayedMatrixStats::rowSums2(getCoverage(bs, type="Cov") >= coverage) >= length(sample.idx))
bs.filtered <- bs[loci.idx, sample.idx]
bs.filtered
head(getCoverage(bs.filtered, type="Cov"))

cat("\n[DM.R] Saving Rdata \t\t\t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
bismark_env <- ls(all=TRUE)
save(list = bismark_env, file = "bismark.RData") 
#load("bismark.RData")

# Distribtuion plots ------------------------------------------------------

cat("\n[DM.R] Plotting Empirical Distribution of CpGs \t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")

# Raw
#plotEmpiricalDistribution(bs, testCovariate=testCovariate)
#plotEmpiricalDistribution(bs, testCovariate=testCovariate, type="Cov", bySample=TRUE)

# Filtered
pdf("Methylation_Distribution.pdf", height = 7.50, width = 11.50)
plotEmpiricalDistribution(bs.filtered, testCovariate=testCovariate)
dev.off()

pdf("Methylation_Coverage.pdf", height = 7.50, width = 11.50)
plotEmpiricalDistribution(bs.filtered, testCovariate=testCovariate, type="Cov", bySample=TRUE)
dev.off()

# DMRs --------------------------------------------------------------------

cat("\n[DM.R] Testing for DMRs with dmrseq \t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")

set.seed(1)
register(MulticoreParam(1))
regions <- dmrseq(bs=bs.filtered,
                  cutoff = 0.05,
                  minNumRegion = 5,
                  maxPerms = 100,
                  testCovariate=testCovariate,
                  adjustCovariate = adjustCovariate,
                  matchCovariate = matchCovariate)

cat("\n[DM.R] Selecting significant DMRs \t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
sum(regions$qval < 0.05)
#sigRegions <- regions[regions$qval < 0.05,]
sigRegions <- regions[regions$pval < 0.05,]
sum(sigRegions$stat > 0) / length(sigRegions)

cat("\n[DM.R] Plotting DMR pie chart \t\t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
pie <- (table(sigRegions$stat < 0))
names(pie) <- c("Hypermethylated", "Hypomethylated")
pdf("HypervsHypo_pie.pdf", height = 8.5, width =11)
pie(pie, labels = c(paste(pie[1], "Hypermethylated", sep = " "), paste(pie[2], "Hypomethylated", sep=" ")), col = c("Red", "Blue"))
dev.off()

cat("\n[DM.R] Extracing raw differnces for DMRs \t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
# Get Raw differences
rawDiff <- meanDiff(bs.filtered, dmrs=regions, testCovariate=testCovariate)
sigRawDiff <- meanDiff(bs.filtered, dmrs=sigRegions, testCovariate=testCovariate)
# Add Raw differences
regions$RawDiff <- rawDiff
sigRegions$RawDiff <- sigRawDiff

cat("\n[DM.R] Exporting DMR and background region information \t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
write.csv(as.data.frame(regions), file="backgroundRegions.csv", row.names = F)
write.csv(as.data.frame(sigRegions), file="DMRs.csv", row.names = F)
write.table(as.data.frame(regions)[1:3], "backgroundRegions.bed", sep ="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(as.data.frame(sigRegions)[1:3], "DMRs.bed", sep ="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

cat("\n[DM.R] Annotating and plotting DMRs \t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
pdf("DMRs.pdf", height = 7.50, width = 11.50)
annoTrack <- getAnnot(genome)
plotDMRs(bs.filtered, regions=sigRegions, testCovariate=testCovariate, annoTrack=annoTrack, qval= F)
dev.off()

cat("\n[DM.R] Saving Rdata \t\t\t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
DMRs_env <- ls(all=TRUE)[!(ls(all=TRUE) %in% bismark_env)]
save(list = DMRs_env, file = "DMRs.RData")
#load("DMRs.RData")

# Individual smoothed values ----------------------------------------------

cat("\n[DM.R] Smoothing individual methylation values \t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
bs.filtered.bsseq <- BSmooth(bs.filtered, mc.cores = cores, verbose = TRUE)

cat("\n[DM.R] Extracting values for DMR heatmap \t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
sig_indiv_smoothed <- data.frame(getMeth(BSseq = bs.filtered.bsseq, regions = sigRegions, type = "smooth", what = "perRegion"))
colnames(sig_indiv_smoothed) <- names
sig_indiv_smoothed_table <- cbind(sigRegions, sig_indiv_smoothed)
write.table(sig_indiv_smoothed_table, "sig_individual_smoothed_DMR_methylation.txt", sep ="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

cat("\n[DM.R] Extracting values for WGCNA \t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
indiv_smoothed <- data.frame(getMeth(BSseq = bs.filtered.bsseq, regions = regions, type = "smooth", what = "perRegion"))
colnames(indiv_smoothed) <- names
indiv_smoothed_table <- cbind(regions, indiv_smoothed)
write.table(indiv_smoothed_table, "background_individual_smoothed_DMR_methylation.txt", sep ="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

cat("\n[DM.R] Creating 20kb windows \t\t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
chrSizes <- seqlengths(goi)
windows <- tileGenome(chrSizes, tilewidth=2e4, cut.last.tile.in.chrom=T)
windows <- keepStandardChromosomes(windows, pruning.mode = "coarse")
windows <- dropSeqlevels(windows, "chrM", pruning.mode = "coarse")
length(seqlevels(windows))
windows

cat("\n[DM.R] Extracting values for 20kb windows \t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
windows_smoothed <- data.frame(getMeth(BSseq = bs.filtered.bsseq, regions = windows, type = "smooth", what = "perRegion"))
colnames(windows_smoothed) <- names
windows_smoothed_table <- cbind(windows, windows_smoothed)
write.table(windows_smoothed_table, "20kb_smoothed_windows.txt", sep ="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

if(genome == "hg38" | genome == "mm10" | genome == "rn6"){
  cat("\n[DM.R] Creating CGi windows \t\t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
  annots = paste(genome,"_cpg_islands", sep="")
  CGi <- build_annotations(genome = genome, annotations = annots)
  CGi <- keepStandardChromosomes(CGi, pruning.mode = "coarse")
  CGi <- dropSeqlevels(CGi, "chrM", pruning.mode = "coarse")
  length(seqlevels(CGi))
  CGi
  
  cat("\n[DM.R] Extracting values for CGi windows \t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
  CGi_smoothed <- data.frame(getMeth(BSseq = bs.filtered.bsseq, regions = CGi, type = "smooth", what = "perRegion"))
  colnames(CGi_smoothed) <- names
  CGi_smoothed_table <- cbind(CGi, CGi_smoothed)
  write.table(CGi_smoothed_table, "CGi_smoothed_windows.txt", sep ="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}

cat("\n[DM.R] Saving Rdata \t\t\t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
bsseq_env <- ls(all=TRUE)[!(ls(all=TRUE) %in% bismark_env) & !(ls(all=TRUE) %in% DMRs_env)]
save(list = bsseq_env, file = "bsseq.RData") 
#load("bsseq.RData")

# PCA of 20 KB windows with CGi -------------------------------------------

cat("\n[DM.R] Tidying 20 kb window data \t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n") 
meth_reorder <- na.omit(windows_smoothed_table[,c(6:length(windows_smoothed_table))])
data <- t(as.matrix(meth_reorder))
# Fix for columns of no variance
#data2 <- as.data.frame(data)
#data3 <- data2[,apply(data2, 2, var, na.rm=TRUE) != 0]
# Titles
stopifnot(sampleNames(bs.filtered.bsseq) == colnames(meth_reorder))
group <- as.character(pData(bs.filtered.bsseq)$Diagnosis)

cat("\n[DM.R] 20 kb window PCA \t\t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
data.pca <- prcomp(data, center = TRUE, scale. = TRUE) 
plot(data.pca, type = "l")
summary(data.pca)

cat("\n[DM.R] Plotting 20 kb window PCA \t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
g <- ggbiplot(data.pca, obs.scale = 1, var.scale = 1, groups = group, ellipse = TRUE, circle = FALSE, var.axes = FALSE, choices = 1:2)
g <- g + scale_color_discrete(name = '')
g <- g + theme_bw(base_size = 25) + geom_point(aes(colour = group), size = 4)
# Change legend position
g <- g + theme(legend.direction = 'vertical',
               legend.position = c(0.125, 0.1),
               legend.text = element_text(size = 12),
               legend.title = element_text(size = 18),
               panel.grid.major = element_blank(), 
               panel.border = element_rect(color = "black", size = 1.25),
               axis.ticks = element_line(size = 1.25), 
               legend.key = element_blank(),
               panel.grid.minor = element_blank()) + guides(col=guide_legend(ncol=2))
#Change title
g <- g + ggtitle("20Kb CpG Windows With CpG Islands")
# Save
print(g)
ggsave("PCA_20kbWindows_withCGi.pdf", plot = last_plot(), device = NULL)

# PCA of CGi windows ------------------------------------------------------

if(genome == "hg38" | genome == "mm10" | genome == "rn6"){
  cat("\n[DM.R] Tidying CGi window data \t\t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
  meth_reorder <- na.omit(CGi_smoothed_table [,c(11:length(CGi_smoothed_table))])
  data <- t(as.matrix(meth_reorder))
  # Fix for columns of no variance
  #data2 <- as.data.frame(data)
  #data3 <- data2[,apply(data2, 2, var, na.rm=TRUE) != 0]
  # Titles
  stopifnot(sampleNames(bs.filtered.bsseq) == colnames(meth_reorder))
  group <- as.character(pData(bs.filtered.bsseq)$Diagnosis)
  
  cat("\n[DM.R] CGi window PCA \t\t\t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n") 
  data.pca <- prcomp(data, center = TRUE, scale. = TRUE) 
  plot(data.pca, type = "l")
  summary(data.pca)
  
  cat("\n[DM.R] Plotting CGi window PCA \t\t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
  g <- ggbiplot(data.pca, obs.scale = 1, var.scale = 1, groups = group, ellipse = TRUE, circle = FALSE, var.axes = FALSE, choices = 1:2)
  g <- g + scale_color_discrete(name = '')
  g <- g + theme_bw(base_size = 25) + geom_point(aes(colour = group), size = 4)
  # Change legend position
  g <- g + theme(legend.direction = 'vertical',
                 legend.position = c(0.125, 0.1),
                 legend.text = element_text(size = 12),
                 legend.title = element_text(size = 18),
                 panel.grid.major = element_blank(), 
                 panel.border = element_rect(color = "black", size = 1.25),
                 axis.ticks = element_line(size = 1.25), 
                 legend.key = element_blank(),
                 panel.grid.minor = element_blank()) + guides(col=guide_legend(ncol=2))
  #Change title
  g <- g + ggtitle("CpG Island Windows")
  # Save
  print(g)
  ggsave("PCA_CGi.pdf", plot = last_plot(), device = NULL)
}

# Heatmap -----------------------------------------------------------------

cat("\n[DM.R] Tidying for heatmap of HCA \t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
# Load smoothed values
matrix <- as.matrix(sig_indiv_smoothed)
# Convert to Percent
matrix <- matrix[,]*100
# Subtract the mean methylation for each row/DMR
data <- sweep(matrix, 1, rowMeans(matrix)) 
# Tidy
data <- as.matrix(data)
colnames(data) <- as.character(pData(bs.filtered.bsseq)$Diagnosis)

cat("\n[DM.R] Plotting heatmap of HCA \t\t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
pdf("heatmap.pdf", height = 8.5, width =11)
heatmap.2(data,
          Rowv= as.dendrogram(hclust(dist(data))),
          scale = c("row"),
          Colv=T,
          col = rev(brewer.pal(11,name="RdBu")),
          margins =c(10,10),
          trace = "none",
          main = paste(nrow(sig_indiv_smoothed),"Differentially Methylated Regions", sep = " "),
          labRow = NA,
          srtCol = 60,
          keysize=0.85,
          key.par = list(cex=0.5),
          key.xlab= "Z-score(% mCG/CG - mean)",
          key.ylab = "Frequency",
          key.title = ""
)
dev.off()

# Prepare files for enrichment analyses -----------------------------------

if(genome == "hg38" | genome == "mm10" | genome == "rn6"){
  cat("\n[CpG_Me] Preparing DMRs files for annotations \t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
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
  
  cat("\n[CpG_Me] Preparing DMRs for HOMER \t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
  dir.create("HOMER")
  write.table((external[,c(1:3)])[external$direction == "Hypermethylated",], "HOMER/DMRs_hyper.bed", quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
  write.table((external[,c(1:3)])[external$direction == "Hypomethylated",], "HOMER/DMRs_hypo.bed", quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
  
  cat("\n[CpG_Me] Preparing background regions for annotations \t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
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
  
  cat("\n[CpG_Me] Building CpG annotations \t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
  annots <- paste(genome,"_cpgs", sep="")
  annotations <- build_annotations(genome = genome, annotations = annots)
  annotations <- keepStandardChromosomes(annotations, pruning.mode = "coarse")
  annotations <- dropSeqlevels(annotations, "chrM", pruning.mode = "coarse")
  length(seqlevels(annotations))
  
  cat("\n[CpG_Me] Annotating DMRs \t\t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
  dm_annotated_CpG <- annotate_regions(
    regions = external,
    annotations = annotations,
    ignore.strand = TRUE,
    quiet = FALSE)
  
  cat("\n[CpG_Me] Annotating background regions \t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
  background_annotated_CpG <- annotate_regions(
    regions = external_bg,
    annotations = annotations,
    ignore.strand = TRUE,
    quiet = FALSE)
  
  cat("\n[CpG_Me] Saving files for GAT \t\t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
  CpGs <- as.data.frame(annotations)
  CpGs <- CpGs[!grepl("_", CpGs$seqnames) ,]
  table(CpGs$seqnames)
  write.table(CpGs[, c(1:3,10)], paste("GAT/",genome,"CpG.bed",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
  
  cat("\n[CpG_Me] Preparing CpG annotation plot \t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
  
  x_order <- c('Hypermethylated','Hypomethylated')
  
  fill_order <- c(
    paste(genome,"_cpg_islands",sep=""),
    paste(genome,"_cpg_shores",sep=""),
    paste(genome,"_cpg_shelves",sep=""),
    paste(genome,"_cpg_inter",sep=""))
  
  CpG_bar <- plot_categorical(
    annotated_regions = dm_annotated_CpG, annotated_random = background_annotated_CpG,
    x='direction', fill='annot.type',
    x_order = x_order, fill_order = fill_order, position='fill',
    plot_title = '',
    legend_title = 'Annotations',
    x_label = '',
    y_label = 'Proportion')
  
  CpG_bar <- CpG_bar + scale_x_discrete(labels=c("All", "Hypermethylated", "Hypomethylated", "Background")) + theme_classic() + 
    theme(axis.text = element_text(size = 25), axis.title = element_text(size = 25),
          strip.text = element_text(size = 25), legend.position = "none", axis.text.x = element_text(angle = 45, hjust=1)) + scale_y_continuous(expand=c(0,0)) 
  
  ggsave("CpG_annotations.pdf", plot = CpG_bar, device = NULL, width = 8.5, height = 11)
  
  # Gene Annotations --------------------------------------------------------
  
  cat("\n[CpG_Me] Building gene region annotations \t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
  annots <- c(paste(genome,"_basicgenes", sep=""),
             paste(genome,"_genes_intergenic", sep=""),
             paste(genome,"_genes_intronexonboundaries", sep=""),
             if(genome == "hg38" | genome == "mm10"){paste(genome,"_enhancers_fantom", sep="")})
  annotations <- build_annotations(genome = genome, annotations = annots)
  annotations <- keepStandardChromosomes(annotations, pruning.mode = "coarse")
  annotations <- dropSeqlevels(annotations, "chrM", pruning.mode = "coarse")
  length(seqlevels(annotations))
  
  cat("\n[CpG_Me] Saving files for GAT \t\t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
  annoFile <- as.data.frame(annotations)
  annoFile <- annoFile[!grepl("_", annoFile$seqnames) ,]
  table(annoFile$seqnames)
  annoFile <- annoFile[, c(1:3,10)]
  
  if(genome == "hg38" | genome == "mm10"){enhancers <-annoFile[annoFile$type == paste(genome,"_enhancers_fantom",sep=""),]}
  promoters <- annoFile[annoFile$type == paste(genome,"_genes_promoters",sep=""),]
  introns <- annoFile[annoFile$type == paste(genome,"_genes_introns",sep=""),]
  boundaries <- annoFile[annoFile$type == paste(genome,"_genes_intronexonboundaries",sep=""),]
  intergenic <- annoFile[annoFile$type == paste(genome,"_genes_intergenic",sep=""),]
  exons <- annoFile[annoFile$type == paste(genome,"_genes_exons",sep=""),]
  fiveUTR <- annoFile[annoFile$type == paste(genome,"_genes_5UTRs",sep=""),]
  threeUTR <- annoFile[annoFile$type == paste(genome,"_genes_3UTRs",sep=""),]
  onetofivekb <- annoFile[annoFile$type == paste(genome,"_genes_1to5kb",sep=""),]
  
  if(genome == "hg38" | genome == "mm10"){write.table(enhancers, "GAT/enhancers.bed", quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")}
  write.table(promoters, "GAT/promoters.bed", quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
  write.table(introns, "GAT/introns.bed", quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
  write.table(boundaries, "GAT/boundaries.bed", quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
  write.table(intergenic, "GAT/intergenic.bed", quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
  write.table(exons, "GAT/exons.bed", quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
  write.table(fiveUTR, "GAT/fiveUTRs.bed", quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
  write.table(threeUTR, "GAT/threeUTRs.bed", quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
  write.table(onetofivekb, "GAT/onetofivekb.bed", quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
  
  cat("\n[CpG_Me] Annotating DMRs \t\t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
  dm_annotated <- annotate_regions(
    regions = external,
    annotations = annotations,
    ignore.strand = TRUE,
    quiet = FALSE)
  
  cat("\n[CpG_Me] Annotating background regions \t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
  background_annotated <- annotate_regions(
    regions = external_bg,
    annotations = annotations,
    ignore.strand = TRUE,
    quiet = FALSE)
  
  cat("\n[CpG_Me] Preparing CpG annotation plot \t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
  
  x_order <- c('Hypermethylated','Hypomethylated')
  
  fill_order <- c(
    if(genome == "hg38" | genome == "mm10"){paste(genome,"_enhancers_fantom",sep="")},
    paste(genome,"_genes_1to5kb",sep=""),
    paste(genome,"_genes_promoters",sep=""),
    paste(genome,"_genes_5UTRs",sep=""),
    paste(genome,"_genes_exons",sep=""),
    paste(genome,"_genes_intronexonboundaries",sep=""),
    paste(genome,"_genes_introns",sep=""),
    paste(genome,"_genes_3UTRs",sep=""),
    paste(genome,"_genes_intergenic",sep=""))
  
  gene_bar <- plot_categorical(
    annotated_regions = dm_annotated, annotated_random = background_annotated,
    x='direction', fill='annot.type',
    x_order = x_order, fill_order = fill_order, position='fill',
    plot_title = '',
    legend_title = 'Annotations',
    x_label = '',
    y_label = 'Proportion')
  
  gene_bar <- gene_bar + scale_x_discrete(labels=c("All", "Hypermethylated", "Hypomethylated", "Background")) + theme_classic() +
    theme(axis.text = element_text(size = 25), axis.title = element_text(size = 25),
          strip.text = element_text(size = 25), legend.position = "none", axis.text.x = element_text(angle = 45, hjust=1)) + scale_y_continuous(expand=c(0,0))
  
  ggsave("generegion_annotations.pdf", plot = region_bar, device = NULL, width = 8.5, height = 11)
}

# GREAT -------------------------------------------------------------------

if(genome=="hg38"){
  cat("\n[DM.R] Obtaining liftOver information for GREAT \t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
  # https://master.bioconductor.org/packages/release/workflows/vignettes/liftOver/inst/doc/liftov.html
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
  
  cat("\n[DM.R] Submitting to GREAT\ \t\t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
  #https://bioconductor.org/packages/release/bioc/vignettes/rGREAT/inst/doc/rGREAT.html
  
  job <- submitGreatJob(sigRegions_liftOver,
                        bg = regions_liftOver,
                        species = "hg19")
  job
  #availableCategories(job)
  #availableOntologies(job)
  tb <- getEnrichmentTables(job, category = c("GO", "Pathway Data"))
  
  cat("\n[DM.R] Saving GREAT enrichment results \t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
  write.xlsx(tb, file = "GREAT_results.xlsx", sep="")
  
  cat("\n[DM.R] Plotting GREAT results \t\t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
  pdf("GREAT.pdf", height = 8.5, width =11)
  par(mfrow = c(1, 3))
  res <- plotRegionGeneAssociationGraphs(job)
  dev.off()
  
  cat("\n[DM.R] Saving GREAT annotations \t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
  res
  write.csv(as.data.frame(res), file="GREATannotations.csv", row.names = F)
}

# Chipseeker --------------------------------------------------------------

# https://bioconductor.org/packages/release/bioc/vignettes/ChIPseeker/inst/doc/ChIPseeker.html

cat("\n[DM.R] Annotating DMRs with ChIPseeker \t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
peakAnno <- annotatePeak(sigRegions,
                         TxDb = TxDb,
                         annoDb = annoDb,
                         overlap = "all")

backgroundAnno <- annotatePeak(regions,
                               TxDb = TxDb,
                               annoDb = annoDb,
                               overlap = "all")

cat("\n[DM.R] ChIPseeker plots \t\t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
# DMR Plots
#pdf("DMRcoverage.pdf", height = 8.5, width =11)
#covplot(sigRegions)
#dev.off()

pdf("Upset.pdf", onefile=FALSE, height = 8.5, width =11)
upsetplot(peakAnno, vennpie=TRUE)
dev.off()

cat("\n[DM.R] Saving ChIPseeker results \t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
annotations <- as.data.frame(peakAnno)
annotations[2] <- annotations[2] - 1
annotations <- annotations[-c(5,17:19)]
write.xlsx(annotations, file = "DMRs_annotated.xlsx", sep= "")

background_annotations <- as.data.frame(backgroundAnno)
background_annotations[2] <- background_annotations[2] - 1
background_annotations <- background_annotations[-c(5,17:19)]
write.xlsx(background_annotations, file = "background_annotated.xlsx", sep= "")

# CMplot ------------------------------------------------------------------

#https://github.com/YinLiLin/R-CMplot

cat("\n[DM.R] Tidying for Manhattan and QQ plots \t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
Manhattan <- as.data.frame(sort(as.GRanges(backgroundAnno), ignore.strand=TRUE))[c("SYMBOL","seqnames", "start", "pval")]
Manhattan$seqnames <- substring(Manhattan$seqnames, 4)

cat("\n[DM.R] Preparing Colors \t\t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
cols = gg_color_hue(2)

cat("\n[DM.R] Generating Manhattan and QQ plots \t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
CMplot(Manhattan,
       col=cols, 
       plot.type=c("m","q"),
       LOG10=TRUE,
       ylim=NULL,
       threshold=0.05,#c(1e-6,1e-4),
       threshold.lty=c(1,2),
       threshold.lwd=c(1,1),
       threshold.col=c("black","grey"),
       cex = 0.5,
       cex.axis = 0.7,
       amplify=FALSE, 
       chr.den.col=brewer.pal(9, "YlOrRd"),
       bin.size=1e6,
       bin.max=100,
       signal.col=c("red","green"),
       signal.cex=c(1,1),
       signal.pch=c(19,19),
       file="pdf",
       memo="")

# Enrichr -----------------------------------------------------------------

# https://cran.r-project.org/web/packages/enrichR/vignettes/enrichR.html
# Check available databases
#dbs <- listEnrichrDbs()

cat("\n[DM.R] Selecting Enrichr databases \t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
dbs <- c("GO_Biological_Process_2018",
         "GO_Cellular_Component_2018",
         "GO_Molecular_Function_2018",
         "KEGG_2016",
         "Panther_2016",
         "Reactome_2016",
         "RNA-Seq_Disease_Gene_and_Drug_Signatures_from_GEO")

cat("\n[DM.R] Running Enrichr \t\t\t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
GO <- enrichr(annotations$SYMBOL, dbs)
write.xlsx(GO, file = "enrichr.xlsx", sep="")

cat("\n[DM.R] Saving RData \t\t\t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
GO_env <- ls(all=TRUE)[!(ls(all=TRUE) %in% bismark_env) & !(ls(all=TRUE) %in% DMRs_env) & !(ls(all=TRUE) %in% bsseq_env)]
save(list = GO_env, file = "GO.RData") 
#load("GO.RData")

# Blocks ------------------------------------------------------------------

cat("\n[DM.R] Testing for large blocks (PMDs/HMDs) \t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")

register(MulticoreParam(1))
blocks <- dmrseq(bs=bs.filtered,
                 cutoff = 0.05,
                 maxPerms = 100,
                 testCovariate=testCovariate,
                 adjustCovariate = adjustCovariate,
                 matchCovariate = matchCovariate,
                 block = TRUE,
                 minInSpan = 500,
                 bpSpan = 5e4,
                 maxGapSmooth = 1e6,
                 maxGap = 5e3)

cat("\n[DM.R] Selecting significant blocks \t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
#sigBlocks <- blocks[blocks$qval < 0.05,]
sigBlocks <- blocks[blocks$pval < 0.05,]

cat("\n[DM.R] Exporting block and background information \t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
write.csv(as.data.frame(blocks), file="backgroundBlocks.csv", row.names = F)
write.csv(as.data.frame(sigBlocks), file="blocks.csv", row.names = F)
write.table(as.data.frame(blocks)[1:3], "backgroundBlocks.bed", sep ="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(as.data.frame(sigRegions)[1:3], "blocks.bed", sep ="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

cat("\n[DM.R] Annotating and plotting blocks \t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
pdf("Blocks.pdf", height = 7.50, width = 11.50)
annoTrack <- getAnnot(genome)
plotDMRs(bs.filtered, regions=sigBlocks, testCovariate=testCovariate, annoTrack=annoTrack, qval= F)
dev.off()

cat("\n[DM.R] Saving RData \t\t\t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
blocks_env <- ls(all=TRUE)[!(ls(all=TRUE) %in% bismark_env) & !(ls(all=TRUE) %in% DMRs_env) & !(ls(all=TRUE) %in% bsseq_env) & !(ls(all=TRUE) %in% GO_env)]
save(list = blocks_env, file = "Blocks.RData") 
#load("Blocks.RData")

# End ---------------------------------------------------------------------

cat("\n[DM.R] Finishing \t\t\t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
sessionInfo()
rm(list=ls())
cat("\n[DM.R] Done \t\t\t\t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
