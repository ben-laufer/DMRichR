#' globalStats
#' @description Computes the average smoothed global and chromosomal CpG methylation values
#'  for each sample and tests for differences between groups while adjusting for the provided covariates. 
#'  Global methylation differences are tested for using an ANOVA through the \code{\link[stats]{aov}} function.
#'  The chromosomal methylation differences are tested using pairwise comparisons calculated from
#'  contrasts of the factor of interest via the \code{\link[lsmeans]{lsmeans}} package.
#' @param bsseq Smoothed bsseq object with design matrix in pData
#' @param testCovar The factor to test for differences between groups
#' @param adjustCovar The covariate(s) to adjust for between groups
#' @param matchCovar Another covariate to adjust for between groups (for dmrseq compatibility)
#' @return A list of tibbles with smoothed global and chromosomal methylation statsitics
#'  and the values used for the tests
#' @references \url{https://cran.r-project.org/web/packages/broom/vignettes/broom_and_dplyr.html}
#' @references \url{https://www.jstatsoft.org/article/view/v069i01/v69i01.pdf}
#' @import bsseq
#' @import GenomeInfoDb
#' @import GenomicRanges
#' @import tidyverse
#' @import lsmeans
#' @import broom
#' @export globalStats
globalStats <- function(bsseq = bs.filtered.bsseq,
                        testCovariate = testCovariate,
                        adjustCovariate = adjustCovariate,
                        matchCovariate = matchCovariate){
  cat("\n[DMRichR] Global and chromosomal methylation statistics \t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
  
  # Linear model formulas ---------------------------------------------------
  cat("Selecting model...")
  
  if(is.null(adjustCovariate) &
     (is.null(matchCovariate) | (length(levels(matchCovariate))) <= 1)){
    model <- as.formula(paste0("CpG_Avg ~ ", paste(testCovariate)))
    
  }else if(!is.null(adjustCovariate) &
           (is.null(matchCovariate) | (length(levels(matchCovariate))) <= 1)){
    model <- as.formula(paste0("CpG_Avg ~ ", paste(testCovariate, "+"), paste(adjustCovariate, collapse = " + ")))
    
  }else if(is.null(adjustCovariate) &
           (!is.null(matchCovariate) | !(length(levels(matchCovariate))) <= 1)){
    model <- as.formula(paste0("CpG_Avg ~ ", paste(testCovariate, "+"), paste(matchCovariate)))
    
  }else if(!is.null(adjustCovariate) &
           (!is.null(matchCovariate) | !(length(levels(matchCovariate))) <= 1)){
    model <- as.formula(paste0("CpG_Avg ~ ", paste(testCovariate, "+"), paste(adjustCovariate, collapse = " + "), paste(" + ", matchCovariate)))
  }
  cat("Done", "\n")
  glue::glue("The model for global and chromosomal statistics is {model}")
  
  # Global ------------------------------------------------------------------
  print(glue::glue("Testing for global methylation differences..."))
  global <- data.frame(DelayedMatrixStats::colMeans2(getMeth(BSseq = bsseq, type = "smooth", what = "perBase")))
  global$sample <- sampleNames(bsseq)
  names(global) <- c("CpG_Avg", "sample")
  global <- dplyr::as_tibble(cbind(global, data.frame(pData(bsseq))), rownames = NULL)
  
  globalResults <- global %>%
    aov(model, data = .) %>% 
    broom::tidy %>% 
    list("globalAnova" = .,
         "globalInput" = global)
  
  # Chromosomal -------------------------------------------------------------
  print(glue::glue("Testing for chromosomal methylation differences..."))
  grl <- split(bsseq, seqnames(bsseq))
  globalChr <- matrix(ncol = length((seqlevels(grl))), nrow = 1)
  for(i in seq_along(seqlevels(grl))){
    globalChr[i] <- data.frame(DelayedMatrixStats::colMeans2(getMeth(BSseq = grl[[i]], type = "smooth", what = "perBase")))
    names(globalChr)[i] <- seqlevels(grl)[i]
  }
  globalChr$sample <- sampleNames(bsseq)
  globalChr <- dplylr::as_tibble(cbind(globalChr, data.frame(pData(bsseq))), rownames = NULL)
  
  pairWise <- globalChr %>% 
    tidyr::gather(key = chromosome,
                  value = CpG_Avg,
                  -sample,
                  -one_of(colnames(pData(bsseq)))) %>% 
    tidyr::nest(-chromosome) %>% 
    dplyr::mutate(
      pairWise = purrr::map(data, ~ lm(model, data = .x) %>% 
                              lsmeans::ref.grid(data = .x) %>%
                              lsmeans::lsmeans(as.formula(paste("~", testCovar))) %>%
                              pairs(reverse = TRUE) %>%
                              summary()
      )
    ) %>%
    dplyr::select(chromosome, pairWise) %>%
    tidyr::unnest() %>%
    dplyr::mutate(fdr = p.adjust(p.value, method = 'fdr'))

  globalResults <- list("globalStats" = globalResults$globalAnova,
                        "globalInput" = globalResults$globalInput,
                        "chromosomalStats" = pairWise,
                        "chromosomalInput" = globalChr)
                         
  return(globalResults)
  
  }
