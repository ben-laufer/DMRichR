#' processBismark
#' @description Process bismark cytosine reports into bsseq objects with design matrix pData
#' @param files List of cytosine report file paths (should only contain reports you are working with for this script)
#' @param names Ordered character vector of sample names
#' @param meta Design matrix table with sample name in the Name column 
#' @param groups Factor of interest (testCovariate)
#' @param mc.cores Number of cores to use
#' @import bsseq
#' @export processBismark
processBismark <- function(files = list.files(path=getwd(), pattern="*.txt.gz"),
                           names =  gsub( "_.*$","", list.files(path=getwd(), pattern="*.txt.gz")),
                           meta = read.csv("sample_info.csv", header = TRUE),
                           groups = testCovariate,
                           Cov = coverage,
                           mc.cores = cores){
  cat("\n[DMRichR] Loading Bismark cytosine reports \t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
  bs <- read.bismark(files = files,
                     sampleNames = names,
                     rmZeroCov = TRUE,
                     strandCollapse = TRUE,
                     fileType = "cytosineReport",
                     verbose = TRUE,
                     mc.cores = mc.cores)
  
  message("Assigning sample metadata...")
  meta <- meta[order(match(meta[,1],names)),]
  stopifnot(sampleNames(bs) == meta$Name)
  pData(bs) <- cbind(pData(bs), meta[2:length(meta)])
  print(pData(bs))
  
  message("Filtering CpGs for coverage...")
  message("Before filtering...")
  bs <- keepStandardChromosomes(bs, pruning.mode = "coarse")
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
