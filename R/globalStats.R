#' globalStats
#' @title Test for global methylation differences
#' @description Computes the average smoothed global CpG methylation values
#'  for each sample and tests for differences between groups while adjusting for the provided 
#'  covariates. CpG island testing is performed for the human, mouse, and rat genomes.
#'  Global methylation and CpG island differences are tested for using an ANOVA through the 
#'  \code{\link[stats]{aov}} function. 
#' @param bs.filtered.bsseq Smoothed \code{bsseq} object with design matrix in pData
#' @param genome Character specifying the genome of interest
#' @param testCovariate The factor to test for differences between groups
#' @param adjustCovariate The covariate(s) to adjust for between groups
#' @param matchCovariate Another covariate to adjust for between groups (for dmrseq compatibility)
#' @return A list of tibbles with smoothed global and CpG island methylation statistics
#'  and the values used for the tests
#' @importFrom magrittr %>%
#' @importFrom DelayedMatrixStats colMeans2
#' @importFrom dplyr as_tibble mutate select
#' @importFrom broom tidy
#' @importFrom tidyr nest unnest gather
#' @importFrom purrr map
#' @importFrom bsseq getMeth
#' @importClassesFrom bsseq BSseq 
#' @importMethodsFrom bsseq pData sampleNames seqnames
#' @import GenomicRanges
#' @importFrom GenomeInfoDb keepStandardChromosomes
#' @importFrom glue glue
#' @importFrom stats aov as.formula
#' @export globalStats
#' 
globalStats <- function(bs.filtered.bsseq = bs.filtered.bsseq,
                        genome = genome,
                        testCovariate = testCovariate,
                        adjustCovariate = NULL,
                        matchCovariate = NULL){
  
  cat("\n[DMRichR] Global and CpG island methylation statistics \t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
  
  # Linear model formulas ---------------------------------------------------
  
  cat("Selecting model...")
  
  if(is.null(adjustCovariate) &
     (is.null(matchCovariate) | (length(levels(matchCovariate))) <= 1)){
    model <- as.formula(paste0("CpG_Avg ~ ", paste(testCovariate)))
    
  }else if(!is.null(adjustCovariate) &
           (is.null(matchCovariate) | (length(levels(matchCovariate))) <= 1)){
    model <- as.formula(paste0("CpG_Avg ~ ", paste(testCovariate, "+"),
                               paste(adjustCovariate, collapse = " + ")))
    
  }else if(is.null(adjustCovariate) &
           (!is.null(matchCovariate) | !(length(levels(matchCovariate))) <= 1)){
    model <- as.formula(paste0("CpG_Avg ~ ", paste(testCovariate, "+"), paste(matchCovariate)))
    
  }else if(!is.null(adjustCovariate) &
           (!is.null(matchCovariate) | !(length(levels(matchCovariate))) <= 1)){
    model <- as.formula(paste0("CpG_Avg ~ ", paste(testCovariate, "+"),
                               paste(adjustCovariate, collapse = " + "), paste(" + ", matchCovariate)))
  }
  cat("Done", "\n")
  cat(paste("The model for globalStats is", paste(capture.output(print(model))[1], collapse= ' ')), "\n")
  
  # Global ------------------------------------------------------------------
  
  cat("Testing for global methylation differences...")
  
  global <- data.frame(DelayedMatrixStats::colMeans2(getMeth(BSseq = bs.filtered.bsseq,
                                                             type = "smooth",
                                                             what = "perBase"),
                                                     na.rm = TRUE))
  global$sample <- sampleNames(bs.filtered.bsseq)
  names(global) <- c("CpG_Avg", "sample")
  global <- dplyr::as_tibble(cbind(global, data.frame(pData(bs.filtered.bsseq))), rownames = NULL)
  
  globalResults <- global %>%
    aov(model, data = .) %>% 
    broom::tidy()
  
  cat("Done", "\n")
  
  # CpG island --------------------------------------------------------------
  
  if(genome %in% c("hg38", "hg19", "mm10", "mm9", "rheMac10", "rheMac8", "rn6", "danRer11", "galGal6",
                   "bosTau9", "panTro6", "dm6", "susScr11", "canFam3")){
    
    print(glue::glue("Performing CpG island methylation level testing for {genome}"))
    
    CGi <- genome %>%
      DMRichR::getCpGs() %>% 
      plyranges::filter(type == "islands") %>% 
      bsseq::getMeth(BSseq = bs.filtered.bsseq,
                     regions = .,
                     type = "smooth",
                     what = "perRegion") %>%
      na.omit() %>% 
      DelayedMatrixStats::colMeans2() %>%
      as.data.frame()
    
    CGi$sample <- sampleNames(bs.filtered.bsseq)
    names(CGi) <- c("CpG_Avg", "sample")
    CGi <- dplyr::as_tibble(cbind(CGi, data.frame(pData(bs.filtered.bsseq))), rownames = NULL)
    
    CGiResults <-  CGi %>%
      aov(model, data = .) %>% 
      broom::tidy() 
    
    cat("Done", "\n")
  }
  
  # Return ------------------------------------------------------------------

  if(genome %in% c("hg38", "hg19", "mm10", "mm9", "rheMac10", "rheMac8", "rn6", "danRer11", "galGal6",
                   "bosTau9", "panTro6", "dm6", "susScr11", "canFam3")){
    
    cat("Returning list of global and chromosomal methylation statistics...")
    
    globalResults <- list("globalInput" = global,
                          "globalStats" = globalResults,
                          "CGiInput" = CGi,
                          "CGiStats" = CGiResults)
    
    cat("Done", "\n")
    
  }else{
    
  cat("Returning list of global and chromosomal methylation statistics...")
    
  globalResults <- list("globalInput" = global,
                        "globalStats" = globalResults)
  
  cat("Done", "\n")
  
  }
                         
  return(globalResults)
  
}
