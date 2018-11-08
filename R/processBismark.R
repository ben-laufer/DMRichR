#' processBismark
#' @description Process bismark cytosine reports into bsseq objects with design matrix pData
#' @param files List of cytosine report file paths
#' @param meta Design matrix table with sample name in the Name column 
#' @param groups Factor of interest (testCovariate)
#' @param Cov Coverage cutoff (1x recommended)
#' @param nThread Number of cores to use
#' @import bsseq
#' @export processBismark
processBismark <- function(files = list.files(path=getwd(), pattern="*.txt.gz"),
                           meta = read.csv("sample_info.csv", header = TRUE),
                           groups = testCovariate,
                           Cov = coverage,
                           nThread = cores){
  cat("\n[DMRichR] Loading Bismark cytosine reports \t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
  files.idx <- pmatch(meta$Name, files)
  files <- files[files.idx]
  names <- as.data.frame(gsub( "_.*$","", files[files.idx]))
  colnames(names) <- "Name"
  rownames(names) <- names[,1]
  names[,1] <- NULL
  
  bs <- read.bismark(files = files,
                     colData = names,
                     rmZeroCov = TRUE,
                     strandCollapse = TRUE,
                     verbose = TRUE,
                     nThread = nThread)
  
  message("Assigning sample metadata...")
  meta <- meta[order(match(meta[,1],rownames(names))),]
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
