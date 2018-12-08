#' globalStats
#' @description Computes the average smoothed global and chromosomal CpG methylation values
#'  for each sample and tests for differences between groups while adjusting for the provided covariates. 
#'  Global methylation differences are tested for using an ANOVA through the \code{\link[stats]{aov}} function.
#'  The chromosomal pairwise contrasts for the factor of interest are computed using the
#'   \code{\link[lsmeans]{lsmeans}} package.
#' @param bsseq Smoothed bsseq object with design matrix in pData
#' @param testCovar The factor to test for differences between groups
#' @param adjustCovar The covariate(s) to adjust for between groups
#' @param matchCovar Another covariate to adjust for between groups (for dmrseq compatibility)
#' @return A list of tibbles with smoothed global and chromosomal methylation statsitics
#'  and the values used for the tests
#' @references \url{https://cran.r-project.org/web/packages/broom/vignettes/broom_and_dplyr.html}
#' @import bsseq
#' @import GenomeInfoDb
#' @import GenomicRanges
#' @import tidyverse
#' @import lsmeans
#' @import broom
#' @export globalStats
globalStats <- function(bsseq = bs.filtered.bsseq,
                        testCovar = testCovariate,
                        adjustCovar = adjustCovariate,
                        matchCovar = matchCovariate){
  cat("\n[DMRichR] Global and chromosomal methylation statistics \t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
  
  # Linear model formulas ---------------------------------------------------
  message("Selecting model...")
  
  if(is.null(adjustCovar) &
     (is.null(matchCovar) | (length(levels(matchCovar))) <= 1)){
    model <- as.formula(paste0("CpG_Avg ~ ", paste(testCovar)))
    
  }else if(!is.null(adjustCovar) &
           (is.null(matchCovar) | (length(levels(matchCovar))) <= 1)){
    model <- as.formula(paste0("CpG_Avg ~ ", paste(testCovar, "+"), paste(adjustCovar, collapse = " + ")))
    
  }else if(is.null(adjustCovar) &
           (!is.null(matchCovar) | !(length(levels(matchCovar))) <= 1)){
    model <- as.formula(paste0("CpG_Avg ~ ", paste(testCovar, "+"), paste(matchCovar)))
    
  }else if(!is.null(adjustCovar) &
           (!is.null(matchCovar) | !(length(levels(matchCovar))) <= 1)){
    model <- as.formula(paste0("CpG_Avg ~ ", paste(testCovar, "+"), paste(adjustCovar, collapse = " + "), paste(" + ", matchCovar)))
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
                       lsmeans(as.formula(paste("~", testCovar))) %>%
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
