#' methylLearn
#' @description Performs feature selection on significant DMRs (predictors) based on random forest (RF) and support vector machine (SVM)
#' algorithms to generate two lists of DMRs ranked by order of importance. Then finds and annotates the DMRs that overlap between the 
#' top 1% of DMRs in the two DMR ranking lists.
#' @param bsseq Smoothed bsseq object.
#' @param regions Genomic ranges object.
#' @param testCovariate Factor of interest.
#' @param TxDb TxDb annotation package for genome of interest.
#' @param annoDb Character specifying OrgDb annotation package for species of interest.
#' @return A tibble object containing the variable importance ranks and annotations of the most important DMRs, which overlap between the 
#' top 1% of DMRs in the two DMR ranking lists obtained from feature selection methods using RF and SVM algorithms.
#' @references \url{https://www.analyticsvidhya.com/blog/2016/03/select-important-variables-boruta-package/}
#' @import bsseq
#' @import tidyverse
#' @import ChIPseeker
#' @import Boruta
#' @import sigFeature
#' @import gt
#' @export methylLearn
methylLearn <- function(bsseq = bs.filtered.bsseq, 
                        regions = sigRegions, 
                        testCovariate = testCovariate, 
                        TxDb = NA, 
                        annoDb = NA) {

  cat(glue::glue("[DMRichR] Learning features of DMRs for {testCovariate}", "\t", format(Sys.time(), "%d-%m-%Y %X"), "\n"))

  # Tidy data ---------------------------------------------------------------
  data <- getMeth(BSseq = bsseq, regions = regions, type = "smooth", what = "perRegion") %>% 
    as.matrix() %>% 
    t()
  groups <- bsseq %>% 
    pData() %>% 
    dplyr::as_tibble() %>% 
    dplyr::pull(!!testCovariate)
  
  # Create a column to replace names for variable importance
  colnames(data) <- cbind(as.data.frame(seqnames(regions)), ranges(regions)) %>%
    dplyr::as_tibble() %>%
    dplyr::select(value,start,end) %>%
    tidyr::unite("bed", c("value","start","end"), sep = ".") %>%
    as.matrix()
  
  # Add "groups" column to "data" tibble 
  data <- data %>% 
    dplyr::as_tibble() %>% 
    tibble::add_column(groups = groups, .before = 1)
  
  
  # Random forest variable importance ---------------------------------------
  # using Boruta algorithm in "Boruta" package
  getBorutaRanking <- function() {
    cat("\n", "Training random forest (RF) model...")
    set.seed(5)
    borutaTrain <- Boruta(groups ~ ., data = data, doTrace = 0)  
    borutaTrainStats <- attStats(borutaTrain)
    
    borutaRanking <- tibble::tibble(DMR = rownames(borutaTrainStats), 
                                    meanImp = borutaTrainStats$meanImp, 
                                    decision = borutaTrainStats$decision) %>% 
      dplyr::arrange(dplyr::desc(meanImp)) %>% 
      tibble::add_column(Rank = 1:nrow(borutaTrainStats), .before = 1)
    
    #plot(borutaRanking$meanImp)
    #plot(borutaRanking$meanImp[1:percent1])
    
    borutaRanking %>%
      dplyr::select(-decision) %>%
      dplyr::rename("DMR [chr.start.end]" = "DMR") %>%
      gt::gt() %>%
      tab_header(
        title = glue::glue("DMR Importance Ranking using Random Forest (RF)"),
        subtitle = glue::glue("meanImp is mean importance value obtained from Boruta algorithm")
      )
    
    cat("Done")
    return(borutaRanking)
  }
  
  
  # Support vector machine variable importance ------------------------------
  # using recursive feature elimination (RFE) algorithm  & t-statistic in "sigFeature" package
  getSigfeatureRanking <- function() {
    cat("\n", "Training support vector machine (SVM) model...")
    dataMatrix <- data %>% 
      dplyr::select(-groups) %>% 
      as.matrix()
    set.seed(5)
    sigfeatureObject <- sigFeature(dataMatrix, data$groups) 
    sigfeatureRanking <- tibble::tibble(Rank = 1:length(sigfeatureObject), 
                                        DMR = colnames(dataMatrix[, sigfeatureObject]))
    sigfeatureRanking %>% 
      dplyr::rename("DMR [chr.start.end]" = "DMR") %>%
      gt::gt() %>%
      tab_header(
        title = glue::glue("DMR Importance Ranking using Support Vector Machine (SVM)"),
        subtitle = glue::glue("obtained from SVM recursive feature elimination (RFE) algorithm and t-statistic")
      )
    
    cat("Done", "\n")
    return(sigfeatureRanking)
  }  
  
  
  rfRanking <- getBorutaRanking()
  svmRanking <- getSigfeatureRanking()
  
  
  # Find overlapping DMRs (predictors) --------------------------------------
  # Find overlap between top 1% predictors in RF and SVM variable importance lists
  numPredictors <- ncol(data) - 1 
  percent1 <- ceiling(.01 * numPredictors) 
  if (percent1 < 10) {
    overlappingTopDmrs <- intersect(rfRanking$DMR[1:10], svmRanking$DMR[1:10]) 
  } else {
    overlappingTopDmrs <- intersect(rfRanking$DMR[1:percent1], svmRanking$DMR[1:percent1]) 
  }
  overlappingTopDmrs_rfRank <- which(rfRanking$DMR %in% overlappingTopDmrs)
  overlappingTopDmrs_svmRank <- which(svmRanking$DMR %in% overlappingTopDmrs)
  
  
  # Annotate DMRs -----------------------------------------------------------
  annotateDmr <- function(dmrList, rfRank, svmRank) {
    annotatedDmrs <- dmrList %>% 
      strsplit(., split = "[.]") %>%
      as.data.frame() %>%
      t() %>%
      magrittr::set_colnames(c("chr", "start", "end")) %>%
      dplyr::as_tibble() %>%
      GenomicRanges::makeGRangesFromDataFrame(ignore.strand = TRUE,
                                              seqnames.field = "chr",
                                              start.field = "start",
                                              end.field = "end") %>%
      ChIPseeker::annotatePeak(TxDb = TxDb,
                               annoDb = annoDb, 
                               overlap = "all") %>%
      dplyr::as_tibble() %>%  
      tibble::add_column(RF.varImp.rank = rfRank, .before = 1) %>%
      tibble::add_column(SVM.varImp.rank = svmRank, .before = 2) %>%
      dplyr::select(-c(strand, geneChr, geneStart, geneEnd, geneStrand, geneId, transcriptId, distanceToTSS, ENSEMBL))
    return(annotatedDmrs)  
  }
  
  annotatedDmrs <- annotateDmr(overlappingTopDmrs, overlappingTopDmrs_rfRank, overlappingTopDmrs_svmRank)
  
  # HTML table of output 
  annotatedDmrs %>% 
    gt::gt() %>%
    tab_header(
      title = glue::glue("Annotations of Overlapping Significant DMRs"),
      subtitle = glue::glue("DMRs overlap between {percent1} (top 1%) out of {numPredictors} DMRs 
                            in two DMR importance ranking lists obtained from feature selection methods 
                            using random forest (RF) and support vector machine (SVM) algorithms")) %>%
    as_raw_html(inline_css = TRUE) %>%
    write("Machine_Learning_Annotated_DMRs.html")

  return(annotatedDmrs)
}



