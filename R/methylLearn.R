#' methylLearn
#' @description Performs feature selection on significant DMRs (predictors) based on random forest (RF) and support vector machine (SVM)
#' algorithms to generate two lists of DMRs ranked by order of importance. Then finds and annotates DMRs that are common among 
#' the top percent (or top 10 or number of predictors if top percent is too low) of DMRs in the two DMR ranking lists.
#' @param bsseq Smoothed bsseq object.
#' @param regions Genomic ranges object.
#' @param testCovariate Factor of interest.
#' @param TxDb TxDb annotation package for genome of interest.
#' @param annoDb Character specifying OrgDb annotation package for species of interest.
#' @param topPercent Positive integer specifying the top percent of DMRs. Default is 1.
#' @param output Either "all" or "one". Default is "all".
#' If "output" is "all", then returned object is a list containing tibbles of: 
#'      1. full RF variable importance ranking, 
#'      2. full SVM variable importance ranking, 
#'      3. annotated DMRs common among the top percent (or top 10 or number of predictors if top percent is too low) of DMRs in the two DMR ranking lists.
#' If "output" is "one", then returned object is a tibble of the annotated common DMRs.
#' @param saveHtmlReport Either TRUE or FALSE. Default is TRUE.
#' If TRUE, an HTML report with the following is generated:
#'      1. Table of annotated top DMRs from RF DMR importance ranking
#'      2. Table of annotated top DMRs from SVM DMR importance ranking
#'      3. Table of annotated common DMRs
#'      4. Heatmap of each sample of each common DMR
#' If FALSE, no HTML report is generated.
#' @return Refer to output argument. Returned object is either a list of tibbles or one tibble.
#' @references \url{https://www.analyticsvidhya.com/blog/2016/03/select-important-variables-boruta-package/}

#' @importFrom dplyr as_tibble select pull arrange desc
#' @importFrom tidyr unite
#' @importFrom tibble tibble add_column
#' @importFrom magrittr %>% set_colnames
#' @import ChIPseeker
#' @import Boruta
#' @import sigFeature
#' @import gt
#' @importFrom pheatmap pheatmap
#' @importFrom glue glue
#' @import R2HTML
#' @importFrom bsseq getMeth
#' @importClassesFrom bsseq BSseq 
#' @importMethodsFrom bsseq pData
#' @export methylLearn
methylLearn <- function(bsseq = bs.filtered.bsseq, 
                        regions = sigRegions, 
                        testCovariate = testCovariate, 
                        TxDb = NA, 
                        annoDb = NA,
                        topPercent = 1,
                        output = "all",
                        saveHtmlReport = TRUE) {

  cat(glue::glue("[DMRichR] Learning features of DMRs for {testCovariate}", "\t", format(Sys.time(), "%d-%m-%Y %X"), "\n"))
  
  # Tidy data ---------------------------------------------------------------
  tidyData <- function() {
    # Get data
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
  

  # Helper function to split DMR to chr, start, end and add to tibble -------
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
  # Using Boruta algorithm in "Boruta" package
  getRfRanking <- function() {
    cat("\n", "Training random forest (RF) model for DMR ranking...")
    set.seed(5)
    borutaTrainObject <- Boruta::Boruta(groups ~ ., data = data %>% dplyr::select(-sampleID), doTrace = 0)  
    borutaTrainStats <- Boruta::attStats(borutaTrainObject)
    
    rfRanking <- tibble::tibble(DMR = rownames(borutaTrainStats), 
                                meanImp = borutaTrainStats$meanImp, 
                                decision = borutaTrainStats$decision) %>% 
      dplyr::arrange(dplyr::desc(meanImp)) %>% 
      tibble::add_column(Rank = 1:nrow(borutaTrainStats), .before = 1) 
    
    rfRanking <- rfRanking %>% splitDmrs()
    cat("Done.")
    return(rfRanking)
  }
  
 
  # Support vector machine variable importance ------------------------------
  # Using recursive feature elimination (RFE) algorithm & t-statistic in "sigFeature" package
  getSvmRanking <- function() {
    cat("\n", "Training support vector machine (SVM) model for DMR ranking...")
    dataMatrix <- data %>% 
      dplyr::select(-c(groups, sampleID)) %>% 
      as.matrix()
    set.seed(5)
    sigfeatTrainObject <- sigFeature::sigFeature(dataMatrix, data$groups) 
    
    svmRanking <- tibble::tibble(Rank = 1:length(sigfeatTrainObject), 
                                 DMR = colnames(dataMatrix[, sigfeatTrainObject]))
    
    svmRanking <- svmRanking %>% splitDmrs()
    cat("Done.")
    return(svmRanking) 
  }  
  

  # Find common DMRs (predictors) --------------------------------------
  # Find DMRs common among top percent (or top 10 or number of predictors if top percent is less than 10) 
  #     of predictors in RF and SVM variable importance lists
  getCommonDmrs <- function() {
    numPredictors <- ncol(data) - 2 
    numTopPercent <- ceiling(topPercent * .01 * numPredictors)
    allCommonFlag <- 0
    
    # Top percent of predictors is less than 10 AND number of predictors is at least 10
    if (numTopPercent < 10 && numPredictors >= 10) {
      commonDmrs <- intersect(rfRanking$DMR[1:10], svmRanking$DMR[1:10]) 
      cat("\n", glue::glue("{topPercent}% of {numPredictors} total DMRs is {numTopPercent} and less than 10."))
      cat("\n", glue::glue("Finding common DMRs in top 10 DMRs (instead of in top {topPercent}%) in RF and SVM predictor ranking lists."))
      case <- 10
      if(length(commonDmrs) == 10) { allCommonFlag <- 1 }
    } 
    # Top percent of predictors is less than 10 AND number of predictors is less than 10
    else if (numTopPercent < 10 && numPredictors < 10) {
      commonDmrs <- intersect(rfRanking$DMR[1:numPredictors], svmRanking$DMR[1:numPredictors]) 
      cat("\n", glue::glue("{topPercent}% of {numPredictors} total DMRs is {numTopPercent} and less than 10."))
      cat("\n", glue::glue("Finding common DMRs in top {numPredictors} DMRs (instead of in top {topPercent}%) in RF and SVM predictor ranking lists."))
      case <- numPredictors
      if(length(commonDmrs) == numPredictors) { allCommonFlag <- 1 }
    } 
    # Top percent of predictors is at least 10
    else {
      cat("\n", glue::glue("{topPercent}% of {numPredictors} total DMRs is {numTopPercent} and at least 10."))
      cat("\n", glue::glue("Finding common DMRs in top {topPercent}% of DMRs in RF and SVM predictor ranking lists."))
      commonDmrs <- intersect(rfRanking$DMR[1:numTopPercent], svmRanking$DMR[1:numTopPercent]) 
      case <- numTopPercent
      if(length(commonDmrs) == numTopPercent) { allCommonFlag <- 1 }
    }
    
    # No common DMRs
    if(length(commonDmrs) == 0) {
      cat("\n", glue::glue("There were 0 common DMRs. Rerun with a higher topPercent value for a greater number of common DMRs."))
    }
    
    # To find correct RF/SVM rank if all selected DMRs are the same in both RF and SVM lists
    if(allCommonFlag == 1) {
      commonDmrsRfRank <- numeric()
      commonDmrsSvmRank <- numeric()
      
      for(i in 1:length(commonDmrs)) {
        # rfRank 
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
    # RF / SVM rank otherwise
    } else {
      commonDmrsRfRank <- which(rfRanking$DMR %in% commonDmrs)
      commonDmrsSvmRank <- which(svmRanking$DMR %in% commonDmrs)
    }
    
    cat("\n")
    return(list(dmrs = commonDmrs, rfRank = commonDmrsRfRank , svmRank = commonDmrsSvmRank, case = case))
  }
 
  
  # Annotate DMRs -----------------------------------------------------------
  annotateDmr <- function(dmrList, rfRank, svmRank, type) {
    cat(" Beginning DMR annotation...", "\n")
    
    annotatedDmrs <- dmrList %>% 
      # Set up DMR list before annotating 
      strsplit(., split = "[.]") %>%
      as.data.frame() %>%
      t() %>%
      magrittr::set_colnames(c("chr", "start", "end")) %>%
      dplyr::as_tibble() %>%
      # Annotate DMRs  
      GenomicRanges::makeGRangesFromDataFrame(ignore.strand = TRUE,
                                              seqnames.field = "chr",
                                              start.field = "start",
                                              end.field = "end") %>%
      ChIPseeker::annotatePeak(TxDb = TxDb,
                               annoDb = annoDb,
                               overlap = "all") %>%
      dplyr::as_tibble() %>%  
      # Add RF rank, SVM rank or both depending on "type"
      {if(type == "rf") tibble::add_column(., rank = rfRank, .before = 1) else .} %>%
      {if(type == "svm") tibble::add_column(., rank = svmRank, .before = 1) else .} %>%
      {if(type == "common") tibble::add_column(., RF.rank = rfRank, .before = 1) else .} %>%
      {if(type == "common") tibble::add_column(., SVM.rank = svmRank, .before = 2) else .} %>%
      # Select only relevant columns
      dplyr::select(-c(strand, geneChr, geneStart, geneEnd, geneStrand, geneId, transcriptId, distanceToTSS, ENSEMBL))
    
    return(annotatedDmrs)  
  }
  

  # Generate HTML table of annotated DMRs -----------------------------------
  getAnnotatedDmrsHtml <- function(annotatedDmrs, type) {
    # Specify titles and subtitles
    if(type == "rf") {
      title <- gt::md("**Annotations of top DMRs from random forest (RF) DMR importance ranking**")
      subtitle <- gt::md("RF DMR importance ranking from Boruta algorithm in *Boruta* package")
    } else if (type == "svm") {
      title <- gt::md("**Annotations of top DMRs from support vector machine (SVM) DMR importance ranking**")
      subtitle <- gt::md("SVM DMR importance ranking from SVM recursive feature elimination (RFE) algorithm & t-statistic in *sigFeature* package")
    } else {
      title <- gt::md("**Annotations of common DMRs**")
      subtitle <- glue::glue("common DMRs from top {commonDmrs$case} DMRs of two DMR importance ranking lists from RF and SVM algorithms")
    }
    
    # Generate HTML of annotated DMRs
    annotatedDmrsHtml <- annotatedDmrs %>% 
      gt::gt() %>%
      gt::tab_header(title = title,
                    subtitle = subtitle)
    
    return(annotatedDmrsHtml)
  }
  

# Generate heatmap of methylation values of each sample of each common DMR --------
  createCommonDmrsHeatmap <- function() {
    # Set up data and labels for heatmap
    heatmapData <- data[, which(colnames(data) %in% commonDmrs$dmrs)] %>% t() 
    colnames(heatmapData) <- data$sampleID 
    annot_col <-  data.frame(testCovariate = data$groups)
    colnames(annot_col) <- testCovariate
    rownames(annot_col) <- colnames(heatmapData)
    
    # Generate heatmap
    commonDmrsHeatmap <- pheatmap::pheatmap(mat = heatmapData, 
                                            angle_col = 45,
                                            border_color = "black", 
                                            main = "Methylation values of each sample of each common DMR",
                                            annotation_col = annot_col,
                                            fontsize = 16,
                                            filename = "./Machine_learning/common_dmrs_heatmap.pdf",
                                            width = 25,
                                            height = 10)
    return(commonDmrsHeatmap)
  }
  
  
  # Create HTML report of:
  #   1. Table of annotated top DMRs from RF DMR importance ranking
  #   2. Table of annotated top DMRs from SVM DMR importance ranking
  #   3. Table of annotated common DMRs
  #   4. Heatmap of each sample of each common DMR
  createHtmlReport <- function() {
    cat(" Generating HTML report...", "\n")
    # annotated top DMRs from RF ranking
    annotatedRfDmrsHtml <- annotateDmr(dmrList = rfRanking$DMR[1:commonDmrs$case], rfRank = rfRanking$Rank[1:commonDmrs$case], svmRank = NULL, type = "rf") %>%
      getAnnotatedDmrsHtml(type = "rf")
    # annotated top DMRs from SVM ranking
    annotatedSvmDmrsHtml <- annotateDmr(dmrList = svmRanking$DMR[1:commonDmrs$case], rfRank = NULL, svmRank = svmRanking$Rank[1:commonDmrs$case], type = "svm") %>%
      getAnnotatedDmrsHtml(type = "svm")
    # annotated common DMRs
    annotatedCommonDmrsHtml <- getAnnotatedDmrsHtml(annotatedCommonDmrs, type = "common")
    # common DMRs heatmap
    createCommonDmrsHeatmap()
    
    
    # output to HTML file
    fileName <- R2HTML::HTMLInitFile(outdir = "./Machine_learning", filename="Machine_learning_report")
    cat("\n<h1 
          align = \"center\"; 
          style= \"font-family: 'Helvetica Neue', Arial, sans-serif; 
          margin-top: 50px;\"> 
          Machine Learning of Significant DMRs
        </h1>", 
        file = fileName)
    R2HTML::HTML("<br>", file = fileName)
    
    R2HTML::HTML(annotatedRfDmrsHtml %>% gt::as_raw_html(inline_css = TRUE), file = fileName)
    R2HTML::HTML("<br>", file = fileName)
    
    R2HTML::HTML(annotatedSvmDmrsHtml %>% gt::as_raw_html(inline_css = TRUE), file = fileName)
    R2HTML::HTML("<br>", file = fileName)
    
    R2HTML::HTML(annotatedCommonDmrsHtml %>% gt::as_raw_html(inline_css = TRUE), file = fileName)
    R2HTML::HTML("<br>", file = fileName)
    
    cat("\n<object data=\"./common_dmrs_heatmap.pdf\" 
            type=\"application/pdf\"
            width=\"1375\"
            height=\"555\">
          alt: <a href=\"./common_dmrs_heatmap.pdf\">
            common_dmrs_heatmap.pdf
          </a>
        </object>",
      file = fileName, append = TRUE)
    R2HTML::HTMLEndFile()
    
    cat(" Done generating HTML report.", "\n")
  }

       
  # Run functions
  data <- tidyData()
  rfRanking <- getRfRanking()
  svmRanking <- getSvmRanking()
  commonDmrs <- getCommonDmrs()
  
  # If no common DMRs, "annotatedCommonDmrs" contains a string message
  if(length(commonDmrs$dmrs) == 0) {
    cat("\n", glue::glue("No annotations due to no common DMRs. Rerun with a higher topPercent value for a greater number of common DMRs."))
    annotatedCommonDmrs <- "No annotations due to no common DMRs." 
  } else {
    annotatedCommonDmrs <- annotateDmr(dmrList = commonDmrs$dmrs, 
                                       rfRank = commonDmrs$rfRank, 
                                       svmRank = commonDmrs$svmRank, 
                                       type = "common")
  }
  
  # If output is "all", "result" is a list of 3 tibbles. If output is "one", "result" is one tibble
  if(output == "all") {
    result <- list("RF ranking" = rfRanking %>% dplyr::select(Rank, DMR, chr, start, end),
                   "SVM ranking" = svmRanking %>% dplyr::select(Rank, DMR, chr, start, end),
                   "Annotated common DMRs" = annotatedCommonDmrs)
  } else {
    result <- annotatedCommonDmrs
  } 
  
  # If "saveHtmlReport" is TRUE, generate HTML report
  if (saveHtmlReport == TRUE) {
    if(!dir.exists("./Machine_learning")) {
      dir.create("./Machine_learning")
    } 
    createHtmlReport()
  } 
  return(result)
}
 