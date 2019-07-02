#' methylLearn
#' @description Performs random forest and suport vector machine feature selection on significant DMRs to identify the most important DMRs
#' @param bsseq Smoothed bsseq object
#' @param regions Genomic ranges object
#' @param groups Ordered test covariate information for each sample
#' @return A tibble object of ranked and annotated DMRs 
#' @references \url{https://cran.rstudio.com/web/packages/randomForestExplainer/vignettes/randomForestExplainer.html}
#' @references \url{https://stats.stackexchange.com/questions/82162/cohens-kappa-in-plain-english}
#' @references \url{https://explained.ai/rf-importance/}
#' @references \url{https://machinelearningmastery.com/feature-selection-with-the-caret-r-package/}
#' @references \url{http://ml-tutorials.kyrcha.info/rfe.html}
#' @import bsseq
#' @import tidyverse
#' @export methylLearn


methylLearn <- function(BSseq = bs.filtered.bsseq, regions = sigRegions) {
  
  # Tidy data ---------------------------------------------------------------
  
  data = getMeth(BSseq = bs.filtered.bsseq, regions = sigRegions, type = "smooth", what = "perRegion") %>% as.matrix() %>% t()
  groups = bs.filtered.bsseq %>% pData() %>% as_tibble() %>% pull(!!testCovariate)
  
  # Create a column to replace names for variable importance
  colnames(data) <- cbind(as.data.frame(seqnames(sigRegions)), ranges(sigRegions)) %>%
    as_tibble() %>%
    dplyr::select(value,start,end) %>%
    tidyr::unite("bed", c("value","start","end"), sep = ".") %>%
    as.matrix()
  
  # Add "groups" column to "data" tibble 
  data <- data %>% as_tibble() %>% add_column(groups = groups, .before = 1)
  
  
  # Random forest variable importance ---------------------------------------
  # using Boruta algorithm in "Boruta" package
  
  getBorutaRanking <- function() {
    set.seed(5)
    borutaTrain <- Boruta(groups ~ ., data = data, doTrace = 2)  
    borutaTrainStats <- attStats(borutaTrain)
    borutaRanking <- tibble(DMR = rownames(borutaTrainStats), meanImp = borutaTrainStats$meanImp, decision = borutaTrainStats$decision) %>% 
      arrange(dplyr::desc(meanImp)) %>% 
      add_column(Ranking = 1:nrow(borutaTrainStats), .before = 1)
    return(borutaRanking)
  }
  
  
  # Support vector machine variable importance ------------------------------
  # using RFE & t-statistic in "sigFeature" package
  
  getSigfeatureRanking <- function() {
    dataMatrix <- data %>% dplyr::select(-groups) %>% as.matrix()
    set.seed(5)
    sigfeatureObject <- sigFeature(dataMatrix, data$groups) 
    sigfeatureRanking <- tibble(Ranking = 1:length(sigfeatureObject), DMR = colnames(dataMatrix[, sigfeatureObject]))
    return(sigfeatureRanking)
  }  
  
  rfRanking <- getBorutaRanking()
  svmRanking <- getSigfeatureRanking()
  
  
  # Find overlapping DMRs (predictors) --------------------------------------
  # Find overlap between top 1% predictors in RF and SVM variable importance lists
  numPredictors = ncol(data) - 1 
  percent1 = ceiling(.01 * numPredictors) 
  
  if (percent1 < 10) {
    overlappingTopDmrs <- intersect(rfRanking$DMR[1:10], svmRanking$DMR[1:10]) 
    overlappingTopDmrs_rfRank <- which(rfRanking$DMR %in% overlappingTopDmrs)
    overlappingTopDmrs_svmRank <- which(svmRanking$DMR %in% overlappingTopDmrs)
  } else {
    overlappingTopDmrs <- intersect(rfRanking$DMR[1:percent1], svmRanking$DMR[1:percent1]) 
    overlappingTopDmrs_rfRank <- which(rfRanking$DMR %in% overlappingTopDmrs)
    overlappingTopDmrs_svmRank <- which(svmRanking$DMR %in% overlappingTopDmrs)
  }
  
  # Annotate DMRs -----------------------------------------------------------
  annotateDmr <- function(dmrList, rfRank, svmRank, filename) {
    dmrSplitted <- strsplit(dmrList, split = "[.]")
    toAnnotate <- data.frame(matrix(unlist(dmrSplitted), nrow = length(dmrSplitted), byrow = T))
    colnames(toAnnotate) <- c("chr", "start", "end")
    
    toAnnotate %>% mutate_if(is.factor, as.character) -> toAnnotate
    
    toAnnotate <- toAnnotate %>%
      GenomicRanges::makeGRangesFromDataFrame(ignore.strand = TRUE,
                                              seqnames.field = "chr",
                                              start.field = "start",
                                              end.field = "end") %>%
      ChIPseeker::annotatePeak(TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene,
                               annoDb = "org.Hs.eg.db",
                               overlap = "all") %>%
      as.data.frame() %>%  
      add_column(RF.varImp.rank = rfRank, .before = 1) %>%
      add_column(SVM.varImp.rank = svmRank, .before = 2) %>%
      select(-c(strand, geneChr, geneStart, geneEnd, geneStrand, geneId, transcriptId, distanceToTSS, ENSEMBL)) 
    
    
    toAnnotate %>% openxlsx::write.xlsx(file = filename)
    return(toAnnotate)  
    
  }
  
  annotateDmr(overlappingTopDmrs, overlappingTopDmrs_rfRank, overlappingTopDmrs_svmRank, "DS_Annotated_Overlapping_Top_DMRs.xlsx")
  
}




