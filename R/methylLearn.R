#' methylLearn
#' @description Performs feature selection on significant DMRs (predictors) based on random forest (RF) and support vector machine (SVM)
#' algorithms to generate two lists of DMRs ranked by order of importance. Then finds and annotates the DMRs that overlap between the 
#' top 1% of DMRs in the two DMR ranking lists.
#' @param bsseq Smoothed bsseq object.
#' @param regions Genomic ranges object.
#' @param testCovariate Factor of interest.
#' @param TxDb TxDb annotation package for genome of interest.
#' @param annoDb Character specifying OrgDb annotation package for species of interest.
#' @param topPercent A positive integer for the top percent of DMRs. 
#' @param output Either "all" or "top percent." The default is "all."
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
                        annoDb = NA,
                        topPercent = 1,
                        output = "all") {

  cat(glue::glue("[DMRichR] Learning features of DMRs for {testCovariate}", "\t", format(Sys.time(), "%d-%m-%Y %X"), "\n"))


  # Error checking ----------------------------------------------------------
  # if(length(topPercent = 1) && output = "multiple top percent") {
  #   cat(glue::glue("Warning: There is only one value indicated for topPercent, but the output option specified is {output}."))
  # }
  
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
    
    borutaTrainObject <- Boruta(groups ~ ., data = data, doTrace = 0)  
    borutaTrainStats <- attStats(borutaTrainObject)
    
    rfRanking <- tibble::tibble(DMR = rownames(borutaTrainStats), 
                                meanImp = borutaTrainStats$meanImp, 
                                decision = borutaTrainStats$decision) %>% 
      dplyr::arrange(dplyr::desc(meanImp)) %>% 
      tibble::add_column(Rank = 1:nrow(borutaTrainStats), .before = 1) 
    
    rfRanking <- rfRanking %>% splitDmrs()
  
    #plot(borutaRanking$meanImp)
    #plot(borutaRanking$meanImp[1:percent1])
    
    cat("Done")
    return(rfRanking)
  }
  

  # Support vector machine variable importance ------------------------------
  # using recursive feature elimination (RFE) algorithm  & t-statistic in "sigFeature" package
  getSvmRanking <- function() {
    cat("\n", "Training support vector machine (SVM) model for DMR ranking...")
    dataMatrix <- data %>% 
      dplyr::select(-groups) %>% 
      as.matrix()
    set.seed(5)
    sigfeatTrainObject <- sigFeature(dataMatrix, data$groups) 
    svmRanking <- tibble::tibble(Rank = 1:length(sigfeatTrainObject), 
                                 DMR = colnames(dataMatrix[, sigfeatTrainObject]))
    
    svmRanking <- svmRanking %>% splitDmrs()
    
    cat("Done", "\n")
    return(svmRanking) 
  }  
  

  
  # Full RF/SVM ranking as HTML table
  getRankingHtml <- function(ranking, type) {
    numPredictors <- ncol(data) - 1 
    numTopPercent <- ceiling(topPercent * .01 * numPredictors)
    if(type == "rf") {
      title <- "Random Forest (RF)"
      subtitle <- "meanImp is mean importance value obtained from Boruta algorithm"
    } else {
      title <- "Support Vector Machine (SVM)"
      subtitle <- "obtained from SVM recursive feature elimination (RFE) algorithm and t-statistic"
    }
    
    # html of full RF/WVM ranking
    htmlRanking <- ranking[1:numTopPercent,] %>%
      dplyr::select(Rank, DMR, chr, start, end) %>%
      dplyr::rename("DMR [chr.start.end]" = "DMR") %>%
      gt::gt() %>%
      tab_header(
        title = glue::glue("DMR Importance Ranking using {title}"),
        subtitle = glue::glue({subtitle})
      )
    return(htmlRanking)
  }
  
  

  
  
  # Find overlapping DMRs (predictors) --------------------------------------
  # Find overlap between top percent of predictors in RF and SVM variable importance lists
  getCommonDmrs <- function() {
    numPredictors <- ncol(data) - 1 
    numTopPercent <- ceiling(topPercent * .01 * numPredictors)
    
    if (numTopPercent < 10 && numPredictors >= 10) {
      overlappingTopDmrs <- intersect(rfRanking$DMR[1:10], svmRanking$DMR[1:10]) 
      cat(glue::glue("{topPercent}% of {numPredictors} total DMRs is {numTopPercent} and less than 10."), "\n")
      cat(glue::glue("Finding common DMRs in top 10 DMRs (instead of in top {topPercent}%) in RF and SVM predictor ranking lists."))
    } else if (numTopPercent < 10 && numPredictors < 10) {
      overlappingTopDmrs <- intersect(rfRanking$DMR[1:numPredictors], svmRanking$DMR[1:numPredictors]) 
      cat(glue::glue("{topPercent}% of {numPredictors} total DMRs is {numTopPercent} and less than 10."), "\n")
      cat(glue::glue("Finding common DMRs in top {numPredictors} DMRs (instead of in top {topPercent}%) in RF and SVM predictor ranking lists."))
    } else {
      cat(glue::glue("Finding common DMRs in top {topPercent}% of DMRs in RF and SVM predictor ranking lists."), "\n")
      overlappingTopDmrs <- intersect(rfRanking$DMR[1:numTopPercent], svmRanking$DMR[1:numTopPercent]) 
      cat(glue::glue("There were {length(overlappingTopDmrs)} common DMRs."))
    }
    
    if(length(overlappingTopDmrs) == 0) {
      cat(glue::glue("There were 0 common DMRs. Rerun with a higher topPercent value for a greater number of common DMRs."))
    }
    overlappingTopDmrs_rfRank <- which(rfRanking$DMR %in% overlappingTopDmrs)
    overlappingTopDmrs_svmRank <- which(svmRanking$DMR %in% overlappingTopDmrs)
    cat("\n")
    
    return(list(dmrs = overlappingTopDmrs, rfRank = overlappingTopDmrs_rfRank, svmRank = overlappingTopDmrs_svmRank))
    
    
  }


  # Heatmap showing correlation of overlapping DMRs -------------------------
  # overlappingTopDmrs
  # head(data)
  # dim(data)
  # data[, which(colnames(data) %in% overlappingTopDmrs)]
  
  
  #pheatmap::pheatmap(data[, which(colnames(data) %in% overlappingTopDmrs)], angle_col = 45)
 
  # Annotate DMRs -----------------------------------------------------------
  annotateDmr <- function(dmrList, rfRank, svmRank) {
    # dmrList <- commonDmrs$dmrs
    # rfRank <- commonDmrs$rfRank
    # svmRank <- commonDmrs$svmRank
    cat("start annotateDmr")
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
      ChIPseeker::annotatePeak(TxDb = TxDb, #TxDb.Hsapiens.UCSC.hg38.knownGene, #TxDb for final
                               annoDb = annoDb, #"org.Hs.eg.db", #annoDb for final
                               overlap = "all") %>%
      dplyr::as_tibble() %>%  
      tibble::add_column(RF.varImp.rank = rfRank, .before = 1) %>%
      tibble::add_column(SVM.varImp.rank = svmRank, .before = 2) %>%
      dplyr::select(-c(strand, geneChr, geneStart, geneEnd, geneStrand, geneId, transcriptId, distanceToTSS, ENSEMBL))
    
    cat("end annotateDmr")
    return(annotatedDmrs)  
  }
  
  # HTML table of output 
  getAnnotatedDmrsHtml <- function(annotatedDmrs) {
    numPredictors <- ncol(data) - 1 
    numTopPercent <- ceiling(topPercent * .01 * numPredictors)
    htmlAnnotatedDmrs <- annotatedDmrs %>% 
      gt::gt() %>%
      gt::tab_header(
        title = glue::glue("Annotations of Overlapping Significant DMRs"),
        subtitle = glue::glue("DMRs overlap between {numTopPercent} (top {topPercent}%) out of {numPredictors} DMRs 
                              in two DMR importance ranking lists obtained from feature selection methods 
                              using random forest (RF) and support vector machine (SVM) algorithms")) 
    return(htmlAnnotatedDmrs)
    # %>%
    # as_raw_html(inline_css = TRUE) %>%
    # write("Machine_Learning_Annotated_DMRs.html")
    
    #annotatedDmrsTable
    
    #(annotatedDmrs$annotation)
  }
  
  getAnnotationCategoryPlot <- function(annotatedDmrs) {
    annotations <- tibble("Annotation" = c("Promoter",
                                           "5' UTR",
                                           "3' UTR",
                                           "Exon",
                                           "Intron",
                                           "Intergenic",
                                           "Downstream"),
                          "Count" = c(str_count(annotatedDmrs$annotation, "Promoter") %>% sum(),
                                      str_count(annotatedDmrs$annotation, "5' UTR") %>% sum(),
                                      str_count(annotatedDmrs$annotation, "3' UTR") %>% sum(),
                                      str_count(annotatedDmrs$annotation, "Exon") %>% sum(),
                                      str_count(annotatedDmrs$annotation, "Intron") %>% sum(),
                                      str_count(annotatedDmrs$annotation, "Intergenic") %>% sum(),
                                      str_count(annotatedDmrs$annotation, "Downstream") %>% sum()) )
    
    dev.new()
    annotationCategoryPlot <- ggplot(annotations, aes(x = Annotation, y = Count)) +
      theme_light() +
      ggtitle("Counts of gene annotation category in overlapping DMRs") +
      xlab("Gene annotation category") +
      geom_bar(stat = "identity", fill = "lightblue") +
      geom_text(aes(label = Count), vjust = 1.4, size = 4)
    
    dev.off()
    return(annotationCategoryPlot)
  }
  
  getHeatmap <- function() {
    dev.new()
    # heatmap shows methylation values of common Dmrs of each sample
    heatmapData <- data[, which(colnames(data) %in% commonDmrs$dmrs)]
    rowName <- character()
    for (i in 1:nrow(heatmapData)) {
      rowName[i] <- glue::glue("Sample {i}")
    }
    heatmapData <- heatmapData %>% 
      as.data.frame() %>%
      magrittr::set_rownames(rowName) 
    
    #mat = data[, which(colnames(data) %in% commonDmrs$dmrs)]
    annotatedDmrsHeatmap <- pheatmap::pheatmap(mat = heatmapData, 
                                               angle_col = 315,
                                               border_color = "black", 
                                               cluster_rows = FALSE, 
                                               cluster_cols = FALSE,
                                               main = "Methylation values of each sample of each common DMR")
    dev.off()
    return(annotatedDmrsHeatmap)
  }

  
  # Run functions
  data <- tidyData() %>% select(c(1:10))
  rfRanking <- getRfRanking()
  svmRanking <- getSvmRanking()
  #rfRankingHtml <- getRankingHtml(ranking = rfRanking, type = "rf")
  #svmRankingHtml <- getRankingHtml(ranking = svmRanking, type = "svm")
  
  commonDmrs <- getCommonDmrs()
  cat("after commonDmrs")
  if(length(commonDmrs$dmr) == 0) {
    annotatedDmrs <- "No annotations due to no common DMRs. Rerun with a higher topPercent value for a greater number of common DMRs."
  } else {
    annotatedDmrs <- annotateDmr(commonDmrs$dmrs, commonDmrs$rfRank, commonDmrs$svmRank)
    #overlappingTopDmrs, overlappingTopDmrs_rfRank, overlappingTopDmrs_svmRank)
  }
  
  #annotatedDmrsHtml <- getAnnotatedDmrsHtml(annotatedDmrs)
  #annotationCategoryPlot <- getAnnotationCategoryPlot(annotatedDmrs)
  #annotatedDmrsHeatmap <- getHeatmap()
  
  if(output == "all") {
    result <- list("RF ranking" = rfRanking %>% select(Rank, DMR, chr, start, end),
                   "SVM ranking" = svmRanking %>% select(Rank, DMR, chr, start, end),
                   "Annotated common DMRS" = annotatedDmrs)
  } else if (output == "top percent") {
    result <- annotatedDmrs
  } 
  
  # end of run functions
  
  
  
  # # HTML report with:
  # rfRankingHtml
  # svmRankingHtml
  # annotatedDmrsHtml
  # annotationCategoryPlot
  # annotatedDmrsHeatmap 
  
  
  return(list(result = result, 
              rfRankingHtml = getRankingHtml(ranking = rfRanking, type = "rf"),
              svmRankingHtml = getRankingHtml(ranking = svmRanking, type = "svm"),
              annotatedDmrsHtml = getAnnotatedDmrsHtml(annotatedDmrs),
              annotationCategoryPlot = getAnnotationCategoryPlot(annotatedDmrs),
              annotatedDmrsHeatmap = getHeatmap()))
}


# 
# 
# library(htmlwidgets)
# library(plotly)
# 
# ggplotly(annotationCounts)
# 
# annotationCounts %>%
#   as_raw_html(inline_css = TRUE) %>%
#   write("report1.html")







# library("DMRichR")
# library("bsseq")
# library("BSgenome.Hsapiens.UCSC.hg38")
# library("TxDb.Hsapiens.UCSC.hg38.knownGene")
# library("org.Hs.eg.db")
# library("tidyverse")
# library("Boruta")
# library("sigFeature")
# library("gt")
# TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene
# annoDb = "org.Hs.eg.db"
# load("../DMRichR_Data/bismark.RData")
# load("../DMRichR_Data/bsseq.RData")
# load("../DMRichR_Data/DMRs.RData")
# load("../DMRichR_Data/GO.RData")
# bsseq = bs.filtered.bsseq
# regions = sigRegions
# testCovariate = testCovariate
# topPercent = 1
# output = "all"


# visualization ideas
# heatmap / how correlated are the top DMRs? doesn't work unless samples are included 
# subset original dataset with the overlapping DMRs to generate correlation heatmap

# color code based on annotation 
# look through DMRichR and find anything to improve on, learn

# TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene
# annoDb = "org.Hs.eg.db"
