#' getGlobal
#' @description Gets smoothed global methylation values
#' @param bsseq Smoothed bsseq object
#' @return Tibble of smoothed global methylation values and phentoype data
#' @import bsseq
#' @import tidyverse
#' @export getGlobal
getGlobal <- function(bsseq = bs.filtered.bsseq){
  if(length(adjustCovariate) == 1){
    cat("\n[DMRichR] Global methylation \t\t\t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
    message("Extracting...")
    global <- data.frame(DelayedMatrixStats::colMeans2(getMeth(BSseq = bsseq, type = "smooth", what = "perBase")))
    global$sample <- sampleNames(bsseq)
    names(global) <- c("CpG_Avg", "sample")
    global <- as.tibble(cbind(global, data.frame(pData(bsseq))), rownames = NULL) %>%
      dplyr::select(sample,
                    CpG_Avg,
                    testCovariate,
                    adjustCovariate,
                    matchCovariate) %>%
      dplyr::rename(sample = "sample",
                    testCovariate = !!testCovariate,
                    adjustCovariate = !!adjustCovariate,
                    matchCovariate = !!matchCovariate)
    return(global)
  }
}

#' getChrom
#' @description Gets smoothed chromosomal methylation values
#' @param bsseq Smoothed bsseq object
#' @return Tibble of smoothed chromosomal methylation values and phentoype data
#' @import bsseq
#' @import GenomeInfoDb
#' @import GenomicRanges
#' @import tidyverse
#' @export getChrom
getChrom <- function(bsseq = bs.filtered.bsseq){
  if(length(adjustCovariate) == 1){
    cat("\n[DMRichR] Chromosomal methylation \t\t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
    message("Extracting...")
    grl <- split(bsseq, seqnames(bsseq))
    global_chr <- matrix(ncol = length((seqlevels(grl))), nrow = 1)
    for(i in seq_along(seqlevels(grl))){
      global_chr[i] <- data.frame(DelayedMatrixStats::colMeans2(getMeth(BSseq = grl[[i]], type = "smooth", what = "perBase")))
      names(global_chr)[i] <- seqlevels(grl)[i]
    }
    global_chr$sample <- sampleNames(bsseq)
    global_chr <- as.tibble(cbind(global_chr, data.frame(pData(bsseq))), rownames = NULL) %>%
      dplyr::rename(sample = "sample",
                    testCovariate = !!testCovariate,
                    adjustCovariate = !!adjustCovariate,
                    matchCovariate = !!matchCovariate) %>%
      dplyr::select(sample,
                    testCovariate,
                    adjustCovariate,
                    matchCovariate,
                    contains("chr")) %>%
      tidyr::gather(key = chromosome,
                    value = CpG_Avg,
                    -sample,
                    -testCovariate,
                    -adjustCovariate,
                    -matchCovariate)
    return(global_chr)
  }
}
