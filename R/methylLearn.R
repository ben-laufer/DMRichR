#' methylLearn
#' @description Performs feature selection on significant DMRs (predictors) based on random forest (RF) and support vector machine (SVM)
#' algorithms to generate two lists of DMRs ranked by order of importance. Then finds and annotates the DMRs that overlap between the 
#' top 1% of DMRs in the two DMR ranking lists.
#' @param bsseq Smoothed bsseq object.
#' @param regions Genomic ranges object.
#' @param testCovariate Factor of interest.
#' @param TxDb TxDb annotation package for genome of interest.
#' @param annoDb Character specifying OrgDb annotation package for species of interest.
#' @param topPercent A positive integer for the top percent of DMRs. Default is 1.
#' @param output Either "all" or "top percent." Default is "all."
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
    cat("\n", "Training random forest (RF) model for DMR ranking...")
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
    cat("\n", "Training support vector machine (SVM) model for DMR ranking...")
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
    
    cat("Done")
    return(sigfeatureRanking)
  }  
  
  
  rfRanking <- getBorutaRanking()
  svmRanking <- getSigfeatureRanking()
  
  # Find overlapping DMRs (predictors) --------------------------------------
  # Find overlap between top percent of predictors in RF and SVM variable importance lists

  numPredictors <- ncol(data) - 1 
  numTopPercent <- ceiling(topPercent * .01 * numPredictors)
  
  if (numTopPercent < 10) {
    overlappingTopDmrs <- intersect(rfRanking$DMR[1:10], svmRanking$DMR[1:10]) 
    cat(glue::glue("{topPercent}% of {numPredictors} total DMRs is {numTopPercent} and less than 10."), "\n")
    cat(glue::glue("Finding common DMRs in top 10 DMRs (not top {topPercent}%) in RF and SVM predictor ranking lists."))
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

  # Heatmap showing correlation of overlapping DMRs -------------------------
  # overlappingTopDmrs
  # head(data)
  # dim(data)
  # data[, which(colnames(data) %in% overlappingTopDmrs)]
  pheatmap::pheatmap(data[, which(colnames(data) %in% overlappingTopDmrs)], angle_col = 45)
  
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
  
  if(length(overlappingTopDmrs) == 0) {
    annotatedDmrs <- "None. No common DMRs."
  } else {
    annotatedDmrs <- annotateDmr(overlappingTopDmrs, overlappingTopDmrs_rfRank, overlappingTopDmrs_svmRank)
  }
  
  if(output == "all") {
    result <- list("RF ranking" = rfRanking[, 1:2],
                   "SVM ranking" = svmRanking,
                   "Annotated common DMRS" = annotatedDmrs)
  } else if (output == "top percent") {
    result <- annotatedDmrs
  } 
  
  # HTML table of output 
  annotatedDmrsTable <- annotatedDmrs %>% 
    gt::gt() %>%
    tab_header(
      title = glue::glue("Annotations of Overlapping Significant DMRs"),
      subtitle = glue::glue("DMRs overlap between {numTopPercent} (top {topPercent}%) out of {numPredictors} DMRs 
                            in two DMR importance ranking lists obtained from feature selection methods 
                            using random forest (RF) and support vector machine (SVM) algorithms")) 
    # %>%
    # as_raw_html(inline_css = TRUE) %>%
    # write("Machine_Learning_Annotated_DMRs.html")
  annotatedDmrsTable
  
  (annotatedDmrs$annotation)
  
  return(result)
}

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
  
annotationCounts <- ggplot(annotations, aes(x = Annotation, y = Count)) + 
  theme_light() + 
  ggtitle("Counts of gene annotation in overlapping DMRs") +
  xlab("Gene annotation") +
  geom_bar(stat = "identity", fill = "lightblue") +
  geom_text(aes(label = Count), vjust = 1.4, size = 4)

annotationCounts

# library("DMRichR")
# library("bsseq")
# library("BSgenome.Hsapiens.UCSC.hg38")
# library("TxDb.Hsapiens.UCSC.hg38.knownGene")
# library("org.Hs.eg.db")
# library("tidyverse")
# library("Boruta")
# library("sigFeature")
# library("gt")
#
# load("../DMRichR_Data/bismark.RData")
# load("../DMRichR_Data/bsseq.RData")
# load("../DMRichR_Data/DMRs.RData")
# load("../DMRichR_Data/GO.RData")

# visualization ideas
# heatmap / how correlated are the top DMRs? doesn't work unless samples are included 
# subset original dataset with the overlapping DMRs to generate correlation heatmap

# color code based on annotation 
# look through DMRichR and find anything to improve on, learn

# TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene
# annoDb = "org.Hs.eg.db"
