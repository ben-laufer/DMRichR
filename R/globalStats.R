#' globalStats
#' @description Obtains smoothed global and chromosomal methylation values and tests for differences between groups while adjusting for the provided covariates
#' @param bsseq Smoothed bsseq object with design matrix in pData
#' @param test The factor to test for differences between groups
#' @param adjust The covariate(s) to adjust for between groups
#' @param match Another covariate to adjust for between groups (for dmrseq compatibility)
#' @return A list of tibbles with smoothed global and chromosomal methylation statsitics and values
#' @references \url{https://cran.r-project.org/web/packages/broom/vignettes/broom_and_dplyr.html}
#' @import bsseq
#' @import GenomeInfoDb
#' @import GenomicRanges
#' @import tidyverse
#' @import lsmeans
#' @import broom
#' @export globalStats
globalStats <- function(bsseq = bs.filtered.bsseq,
                        test = testCovariate,
                        adjust = adjustCovariate,
                        match = matchCovariate){
  cat("\n[DMRichR] Global and chromosomal methylation statistics \t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
  
  # Linear model formulas ---------------------------------------------------
  message("Selecting model...")
  
  if(is.null(adjust) &
     (is.null(match) | (length(levels(match))) <= 1)){
    model <- as.formula(paste0("CpG_Avg ~ ", paste(test)))
    
  }else if(!is.null(adjust) &
           (is.null(match) | (length(levels(match))) <= 1)){
    model <- as.formula(paste0("CpG_Avg ~ ", paste(test, "+"), paste(adjust, collapse = " + ")))
    
  }else if(is.null(adjust) &
           (!is.null(match) | !(length(levels(match))) <= 1)){
    model <- as.formula(paste0("CpG_Avg ~ ", paste(test, "+"), paste(match)))
    
  }else if(!is.null(adjust) &
           (!is.null(match) | !(length(levels(match))) <= 1)){
    model <- as.formula(paste0("CpG_Avg ~ ", paste(test, "+"), paste(adjust, collapse = " + "), paste(" + ", match)))
  }
  
  
  # Global ------------------------------------------------------------------
  message("Testing for global methylation differences...")
  global <- data.frame(DelayedMatrixStats::colMeans2(getMeth(BSseq = bsseq, type = "smooth", what = "perBase")))
  global$sample <- sampleNames(bsseq)
  names(global) <- c("CpG_Avg", "sample")
  global <- as.tibble(cbind(global, data.frame(pData(bsseq))), rownames = NULL)
  
  globalResults <- global %>%
    aov(model, data = .) %>% 
    tidy %>% 
    list("globalAnova" = .,
         "globalInput" = global)
  
  # Chromosomal -------------------------------------------------------------
  message("Testing for chromosomal methylation differences...")
  grl <- split(bsseq, seqnames(bsseq))
  globalChr <- matrix(ncol = length((seqlevels(grl))), nrow = 1)
  for(i in seq_along(seqlevels(grl))){
    globalChr[i] <- data.frame(DelayedMatrixStats::colMeans2(getMeth(BSseq = grl[[i]], type = "smooth", what = "perBase")))
    names(globalChr)[i] <- seqlevels(grl)[i]
  }
  globalChr$sample <- sampleNames(bsseq)
  globalChr <- as.tibble(cbind(globalChr, data.frame(pData(bsseq))), rownames = NULL)
  
  pairWise <- globalChr %>% 
    tidyr::gather(key = chromosome,
                  value = CpG_Avg,
                  -sample,
                  -one_of(colnames(pData(bsseq)))) %>% 
    nest(-chromosome) %>% 
    mutate(
      pairWise = map(data, ~ lm(model, data = .x) %>% 
                       ref.grid(data = .x) %>%
                       lsmeans(as.formula(paste("~", test))) %>%
                       pairs() %>%
                       summary()
                     )
    ) %>%
    dplyr::select(chromosome, pairWise) %>%
    unnest  %>%
    mutate(fdr = p.adjust(p.value, method = 'fdr'))

  globalResults <- list("globalStats" = globalResults$globalAnova,
                        "globalInput" = globalResults$globalInput,
                        "chromosomalStats" = pairWise,
                        "chromosomalInput" = globalChr)
                         
  return(globalResults)
  
  }
