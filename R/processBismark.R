#' processBismark
#' @description Process bismark cytosine reports into bsseq objects with design matrix pData
#' @param files List of cytosine report file paths
#' @param meta Design matrix table with sample name in the Name column 
#' @param groups Factor of interest (testCovariate)
#' @param Cov Coverage cutoff (1x recommended)
#' @param mc.cores Number of cores to use
#' @import bsseq
#' @import openxlsx
#' @export processBismark
processBismark <- function(files = list.files(path = getwd(), pattern = "*.txt.gz"),
                           meta = read.xlsx("sample_info.xlsx", colNames = TRUE),
                           groups = testCovariate,
                           Cov = coverage,
                           mc.cores = cores){
  cat("\n[DMRichR] Loading Bismark cytosine reports \t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
  message("Selecting files...")
  files.idx <- pmatch(meta$Name, files)
  files <- files[files.idx]
  #names <- as.data.frame(gsub( "_.*$","", files[files.idx]))
  #colnames(names) <- "Name"
  #rownames(names) <- names[,1]
  #names[,1] <- NULL
  
  #message("Determining parallelization...")
  #if(mc.cores > 1){
  #  BPPARAM <- BiocParallel::MulticoreParam(workers = floor(mc.cores/4), progressbar = TRUE)
  #  nThread <- floor(mc.cores/floor(mc.cores/4))
  #  message("Parallel processing will be used")
  #}else if(mc.cores == 1){
  #  BPPARAM <- BiocParallel::MulticoreParam(workers = mc.cores, progressbar = TRUE)
  #  nThread <- mc.cores
  #  message("Parallel processing will not be used")
  #}
 
  message("Reading cytosine reports...")
  bs <- read.bismark(files = files,
                     #colData = names,
                     rmZeroCov = FALSE,
                     strandCollapse = TRUE,
                     verbose = TRUE,
                     BPPARAM = MulticoreParam(workers = mc.cores, progressbar = TRUE), # BPPARAM # bpparam()
                     nThread = 1) # nThread
  
  message("Assigning sample metadata...")
  sampleNames(bs) <- gsub( "_.*$","", sampleNames(bs))
  meta <- meta[order(match(meta[,1],sampleNames(bs))),]
  stopifnot(sampleNames(bs) == as.character(meta$Name))
  pData(bs) <- cbind(pData(bs), meta[2:length(meta)])
  print(pData(bs))
  
  message("Filtering CpGs for coverage...")
  message("Before filtering...")
  bs <- GenomeInfoDb::keepStandardChromosomes(bs, pruning.mode = "coarse")
  print(bs)
  print(head(getCoverage(bs, type = "Cov")))
  sample.idx <- which(pData(bs)[[groups]] %in% levels(pData(bs)[[groups]]))
  loci.idx <- which(DelayedMatrixStats::rowSums2(getCoverage(bs, type="Cov") >= Cov) >= length(sample.idx))
  bs.filtered <- bs[loci.idx, sample.idx]
  message("After filtering...")
  print(head(getCoverage(bs.filtered, type = "Cov")))
  print(bs.filtered)
  return(bs.filtered)
}
