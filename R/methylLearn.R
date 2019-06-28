#' methylLearn
#' @description Performs random forest and suport vector machine feature selection on significant DMRs to identify the most important DMRs
#' @param bsseq Smoothed bsseq object
#' @param regions Genomic ranges object
#' @param names Ordered sample names
#' @param groups Ordered test covariate information for each sample
#' @param k.fold number of folds in the cross-validation
#' @param repeats number of repeats in the cross-validation
#' @param ... Additional arguments passed onto randomForest::randomForest()
#' @return A tibble object of ranked and annotated DMRs 
#' @references \url{https://cran.rstudio.com/web/packages/randomForestExplainer/vignettes/randomForestExplainer.html}
#' @references \url{https://stats.stackexchange.com/questions/82162/cohens-kappa-in-plain-english}
#' @references \url{https://explained.ai/rf-importance/}
#' @references \url{https://machinelearningmastery.com/feature-selection-with-the-caret-r-package/}
#' @references \url{http://ml-tutorials.kyrcha.info/rfe.html}
#' @import bsseq
#' @import tidyverse
#' @export methylLearn


# RF variable importance vs. Boruta:
# https://www.listendata.com/2017/05/feature-selection-boruta-package.html
# Boruta:
# https://www.analyticsvidhya.com/blog/2016/03/select-important-variables-boruta-package/
# 
## from source: http://r-statistics.co/Variable-Selection-and-Importance-With-R.html
#
# # Decide if a variable is important or not using Boruta
# boruta_output <- Boruta(ozone_reading ~ ., data=na.omit(inputData), doTrace=2)  # perform Boruta search
# # Confirmed 10 attributes: Humidity, Inversion_base_height, Inversion_temperature, Month, Pressure_gradient and 5 more.
# # Rejected 3 attributes: Day_of_month, Day_of_week, Wind_speed.
# boruta_signif <- names(boruta_output$finalDecision[boruta_output$finalDecision %in% c("Confirmed", "Tentative")])  # collect Confirmed and Tentative variables
# print(boruta_signif)  # significant variables
# #=> [1] "Month"                 "ozone_reading"         "pressure_height"      
# #=> [4] "Humidity"              "Temperature_Sandburg"  "Temperature_ElMonte"  
# #=> [7] "Inversion_base_height" "Pressure_gradient"     "Inversion_temperature"
# #=> [10] "Visibility"
# plot(boruta_output, cex.axis=.7, las=2, xlab="", main="Variable Importance")  # plot variable importance
##
##


library("DMRichR")
library("bsseq")
library("BSgenome.Hsapiens.UCSC.hg38")
library("TxDb.Hsapiens.UCSC.hg38.knownGene")
library("org.Hs.eg.db")
library("tidyverse")
library("Boruta")
library("sigFeature")

load("../DMRichR_Data/bismark.RData")
load("../DMRichR_Data/bsseq.RData")
load("../DMRichR_Data/DMRs.RData")
load("../DMRichR_Data/GO.RData")


methylLearn <- function(data = getMeth(BSseq = bs.filtered.bsseq, regions = sigRegions, type = "smooth", what = "perRegion") %>% as.matrix() %>% t(),
                        groups = bs.filtered.bsseq %>% pData() %>% as_tibble() %>% pull(!!testCovariate)) {

  # Tidy data ---------------------------------------------------------------

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
    # from running boruta.train:
    # There were 34 warnings (use warnings() to see them)
    # 1: In addShadowsAndGetImp(decReg, runs) :
    #    getImp result contains NA(s) or NaN(s); replacing with 0(s), yet this is suspicious.
    borutaTrainStats <- attStats(borutaTrain)
    borutaRanking <- tibble(DMR = rownames(borutaTrainStats), meanImp = borutaTrainStats$meanImp, decision = borutaTrainStats$decision) %>% 
      arrange(dplyr::desc(meanImp)) %>% 
      add_column(Ranking = 1:nrow(borutaTrainStats), .before = 1)
    return(borutaRanking)
    # borutaRankedDmrs <- borutaRanking$DMR
    # return(borutaRankedDmrs)
  }
  

  # Support vector machine variable importance ------------------------------
  # using RFE & t-statistic in "sigFeature" package
  
  getSigfeatureRanking <- function() {
    dataMatrix <- data %>% dplyr::select(-groups) %>% as.matrix()
    sigfeatureObject <- sigFeature(dataMatrix, data$groups) # start 3:21-3:22 
    sigfeatureRanking <- tibble(Ranking = 1:length(sigfeatureObject), DMR = colnames(dataMatrix[, sigfeatureObject]))
     
    return(sigfeatureRanking)
    # sigfeatureRankedDmrs <- colnames(dataMatrix[, sigfeatureRanking]) 
    # return(sigfeatureRankedDmrs)
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
    
    toAnnotate %>%
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
      select(-c(strand, geneChr, geneStart, geneEnd, geneStrand, geneId, transcriptId, distanceToTSS, ENSEMBL)) %>%
      openxlsx::write.xlsx(file = filename)
  }
  
  annotateDmr(overlappingTopDmrs, overlappingTopDmrs_rfRank, overlappingTopDmrs_svmRank, "DS_Annotated_Overlapping_Top_DMRs.xlsx")
 
}
 
# ~ 6 minutes 
methylLearn()










