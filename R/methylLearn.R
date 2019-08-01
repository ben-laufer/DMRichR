#' methylLearn
#' @description Performs feature selection on significant DMRs (predictors) based on random forest (RF) and support vector machine (SVM)
#' algorithms to generate two lists of DMRs ranked by order of importance. Then finds and annotates the common DMRs that overlap between 
#' the top percent (or top 10 or number of predictors if top percent is too low) of DMRs in the two DMR ranking lists.
#' @param bsseq Smoothed bsseq object.
#' @param regions Genomic ranges object.
#' @param testCovariate Factor of interest.
#' @param TxDb TxDb annotation package for genome of interest.
#' @param annoDb Character specifying OrgDb annotation package for species of interest.
#' @param topPercent A positive integer specifying the top percent of DMRs. Default is 1.
#' @param output Either "all" or "top". Default is "all".
#' @return A list containing: result, rfRankingHtml, svmRankingHtml, annotatedDmrsHtml, annotatedDmrsHeatmap.
#' If "output" == "all", then "result" is a list containing tibbles of: 
#'      full variable importance RF ranking, 
#'      full variable importance SVM ranking, 
#'      annotated common DMRs that overlap between the top percent (or top 10 or number of predictors if top percent is too low) of DMRs in the two DMR ranking lists.
#' If "output" == "top", then "result" is a tibble of the annotated common DMRs.
#' "rfRankingHtml" and "svmRankingHtml" are the two DMR ranking lists used to find common DMRs (not the full ranking for all DMRs)
#' "annotatedDmrsHtml" and "annotatedDmrsHeatmap" are an HTML table and heatmap of common DMRs that were annotated.
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
                        annoDb = NA,
                        topPercent = 1,
                        output = "all") {

  cat(glue::glue("[DMRichR] Learning features of DMRs for {testCovariate}", "\t", format(Sys.time(), "%d-%m-%Y %X"), "\n"))
  
  # Tidy data ---------------------------------------------------------------
  tidyData <- function() {
    data <- getMeth(BSseq = bsseq, regions = regions, type = "smooth", what = "perRegion") %>% 
      as.matrix() %>% 
      t()

    # Create a column to replace names for variable importance
    colnames(data) <- cbind(as.data.frame(seqnames(regions)), ranges(regions)) %>%
      dplyr::as_tibble() %>%
      dplyr::select(value,start,end) %>%
      tidyr::unite("bed", c("value","start","end"), sep = ".") %>%
      as.matrix()
    
    # Add "groups" column to "data" tibble 
    data <- data %>% 
      dplyr::as_tibble() %>% 
      tibble::add_column(groups = bsseq %>% 
                           pData() %>% 
                           dplyr::as_tibble() %>% 
                           dplyr::pull(!!testCovariate),
                         .before = 1) %>%
      tibble::add_column(sampleID = bsseq %>%
                           pData() %>%
                           rownames(),
                         .before = 1)

    return(data)
  }
  
  # Helper function to split DMR to chr, start, end and add to tibble
  splitDmrs <- function(ranking) {
    DMRsplit <- ranking$DMR %>%
      strsplit(., split = "[.]") %>%
      as.data.frame() %>%
      t() %>%
      magrittr::set_colnames(c("chr", "start", "end")) %>%
      dplyr::as_tibble()
    
    ranking <- ranking %>%
      tibble::add_column(chr = DMRsplit$chr, .after = 2) %>%
      tibble::add_column(start = DMRsplit$start, .after = 3) %>%
      tibble::add_column(end = DMRsplit$end, .after = 4)
    
    return(ranking)
  }
  
  # Random forest variable importance ---------------------------------------
  # using Boruta algorithm in "Boruta" package
  getRfRanking <- function() {
    cat("\n", "Training random forest (RF) model for DMR ranking...")
    set.seed(5)
    borutaTrainObject <- Boruta(groups ~ ., data = data %>% select(-sampleID), doTrace = 0)  
    borutaTrainStats <- attStats(borutaTrainObject)
    
    rfRanking <- tibble::tibble(DMR = rownames(borutaTrainStats), 
                                meanImp = borutaTrainStats$meanImp, 
                                decision = borutaTrainStats$decision) %>% 
      dplyr::arrange(dplyr::desc(meanImp)) %>% 
      tibble::add_column(Rank = 1:nrow(borutaTrainStats), .before = 1) 
    
    rfRanking <- rfRanking %>% splitDmrs()
    cat("Done")
    return(rfRanking)
  }
  

  # Support vector machine variable importance ------------------------------
  # using recursive feature elimination (RFE) algorithm  & t-statistic in "sigFeature" package
  getSvmRanking <- function() {
    cat("\n", "Training support vector machine (SVM) model for DMR ranking...")
    dataMatrix <- data %>% 
      dplyr::select(-c(groups, sampleID)) %>% 
      as.matrix()
    set.seed(5)
    sigfeatTrainObject <- sigFeature(dataMatrix, data$groups) 
    
    svmRanking <- tibble::tibble(Rank = 1:length(sigfeatTrainObject), 
                                 DMR = colnames(dataMatrix[, sigfeatTrainObject]))
    
    svmRanking <- svmRanking %>% splitDmrs()
    cat("Done")
    return(svmRanking) 
  }  
  

  # Find overlapping/common DMRs (predictors) --------------------------------------
  # Find overlap between top percent of predictors in RF and SVM variable importance lists
  getCommonDmrs <- function() {
    numPredictors <- ncol(data) - 2 
    numTopPercent <- ceiling(topPercent * .01 * numPredictors)
    allCommonFlag <- 0
    
    # top percent of predictors is less than 10 AND number of predictors is at least 10
    if (numTopPercent < 10 && numPredictors >= 10) {
      commonDmrs <- intersect(rfRanking$DMR[1:10], svmRanking$DMR[1:10]) 
      cat("\n", glue::glue("{topPercent}% of {numPredictors} total DMRs is {numTopPercent} and less than 10."))
      cat("\n", glue::glue("Finding common DMRs in top 10 DMRs (instead of in top {topPercent}%) in RF and SVM predictor ranking lists."))
      case <- 10
      if(length(commonDmrs) == 10) { allCommonFlag <- 1 }
    } 
    # top percent of predictors is less than 10 AND number of predictors is less than 10
    else if (numTopPercent < 10 && numPredictors < 10) {
      commonDmrs <- intersect(rfRanking$DMR[1:numPredictors], svmRanking$DMR[1:numPredictors]) 
      cat("\n", glue::glue("{topPercent}% of {numPredictors} total DMRs is {numTopPercent} and less than 10."))
      cat("\n", glue::glue("Finding common DMRs in top {numPredictors} DMRs (instead of in top {topPercent}%) in RF and SVM predictor ranking lists."))
      case <- numPredictors
      if(length(commonDmrs) == numPredictors) { allCommonFlag <- 1 }
    } 
    # top percent of predictors is at least 10
    else {
      cat("\n", glue::glue("{topPercent}% of {numPredictors} total DMRs is {numTopPercent} and at least 10."))
      cat("\n", glue::glue("Finding common DMRs in top {topPercent}% of DMRs in RF and SVM predictor ranking lists."))
      commonDmrs <- intersect(rfRanking$DMR[1:numTopPercent], svmRanking$DMR[1:numTopPercent]) 
      case <- numTopPercent
      if(length(commonDmrs) == numTopPercent) { allCommonFlag <- 1 }
    }
    
    # if no overlapping/common DMRs
    if(length(commonDmrs) == 0) {
      cat("\n", glue::glue("There were 0 common DMRs. Rerun with a higher topPercent value for a greater number of common DMRs."))
    }
    
    # To find correct RF/SVM rank if all selected DMRs are the same in both RF and SVM lists
    if(allCommonFlag == 1) {
      commonDmrsRfRank <- numeric()
      commonDmrsSvmRank <- numeric()
      
      for(i in 1:length(commonDmrs)) {
        # rfRank is already in order from 1:number of common DMRs
        for (rfDmr in rfRanking$DMR) {
          if(commonDmrs[i] == rfDmr) {
            rfRank <- rfRanking$Rank[which(rfRanking$DMR == rfDmr)]
            commonDmrsRfRank <- commonDmrsRfRank %>% append(rfRank)
          }
        }
        # svmRank
        for (svmDmr in svmRanking$DMR) {
          if(commonDmrs[i] == svmDmr) {
            svmRank <- svmRanking$Rank[which(svmRanking$DMR == svmDmr)]
            commonDmrsSvmRank <- commonDmrsSvmRank %>% append(svmRank)
          }
        }
      }
    # RF/SVM rank otherwise
    } else {
      commonDmrsRfRank <- which(rfRanking$DMR %in% commonDmrs)
      commonDmrsSvmRank <- which(svmRanking$DMR %in% commonDmrs)
    }
    
    return(list(dmrs = commonDmrs, rfRank = commonDmrsRfRank , svmRank = commonDmrsSvmRank, case = case))
  }
 
  # Annotate DMRs -----------------------------------------------------------
  annotateDmr <- function(dmrList, rfRank, svmRank) {
    cat("\n", "Beginning DMR annotation...", "\n")
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
  
  # HTML table of output 
  getAnnotatedDmrsHtml <- function(annotatedDmrs) {
    numPredictors <- ncol(data) - 2 
    numTopPercent <- ceiling(topPercent * .01 * numPredictors)
    htmlAnnotatedDmrs <- annotatedDmrs %>% 
      gt::gt() %>%
      gt::tab_header(
        title = glue::glue("Annotations of overlapping significant DMRs"),
        subtitle = glue::glue("DMRs overlap between {numTopPercent} (top {topPercent}%) out of {numPredictors} DMRs 
                              in two DMR importance ranking lists obtained from feature selection methods 
                              using random forest (RF) and support vector machine (SVM) algorithms")) 
    return(htmlAnnotatedDmrs)
  }
  
  
  # RF/SVM ranking of Top percent/predictors as HTML table
  getRankingHtml <- function(ranking, type, case) {
    if(type == "rf") {
      title <- "random forest (RF)"
      subtitle <- "meanImp is mean importance value obtained from Boruta algorithm"
    } else {
      title <- "support vector machine (SVM)"
      subtitle <- "obtained from SVM recursive feature elimination (RFE) algorithm and t-statistic"
    }
    
    # html of full RF/WVM ranking
    htmlRanking <- ranking[1:case,] %>%
      dplyr::select(Rank, DMR, chr, start, end) %>%
      dplyr::rename("DMR [chr.start.end]" = "DMR") %>%
      gt::gt() %>%
      tab_header(
        title = glue::glue("DMR importance ranking using {title}"),
        subtitle = glue::glue({subtitle})
      )
    return(htmlRanking)
  }
  
  getHeatmap <- function() {
    while (dev.cur() > 1) dev.off()
    
    # heatmap shows methylation values of common Dmrs of each sample
    heatmapData <- data[, which(colnames(data) %in% commonDmrs$dmrs)] %>% t() 
    colnames(heatmapData) <- data$sampleID 
    
    annot_col <-  data.frame(testCovariate = data$groups)
    colnames(annot_col) <- testCovariate
    rownames(annot_col) <- colnames(heatmapData)
    
    annotatedDmrsHeatmap <- pheatmap::pheatmap(mat = heatmapData, 
                                               angle_col = 45,
                                               border_color = "black", 
                                               main = "Methylation values of each sample of each common DMR",
                                               annotation_col = annot_col,
                                               fontsize = 14)
    
    annotatedDmrsHeatmap 
    
    return(annotatedDmrsHeatmap)
  }
  
  # Run functions
  data <- tidyData()
  rfRanking <- getRfRanking()
  svmRanking <- getSvmRanking()
  commonDmrs <- getCommonDmrs()
  
  if(length(commonDmrs$dmrs) == 0) {
    cat("\n", glue::glue("No annotations due to no common DMRs. Rerun with a higher topPercent value for a greater number of common DMRs."))
    annotatedDmrs <- "No annotations due to no common DMRs." 
  } else {
    annotatedDmrs <- annotateDmr(commonDmrs$dmrs, commonDmrs$rfRank, commonDmrs$svmRank)
  }
  
  if(output == "all") {
    result <- list("RF ranking" = rfRanking %>% select(Rank, DMR, chr, start, end),
                   "SVM ranking" = svmRanking %>% select(Rank, DMR, chr, start, end),
                   "Annotated common DMRs" = annotatedDmrs)
  } else {
    result <- annotatedDmrs
  } 
  
  # HTML report with:
  #     rfRankingHtml
  #     svmRankingHtml
  #     annotatedDmrsHtml
  #     annotatedDmrsHeatmap 
  
  return(list(result = result, 
              rfRankingHtml = getRankingHtml(ranking = rfRanking, type = "rf", case = commonDmrs$case),
              svmRankingHtml = getRankingHtml(ranking = svmRanking, type = "svm", case = commonDmrs$case),
              annotatedDmrsHtml = getAnnotatedDmrsHtml(annotatedDmrs),
              annotatedDmrsHeatmap = getHeatmap()))
}
