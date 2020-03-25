#' chromHMM
#' @description Perfom enrichment testing against the ChromHMM 15-state model for hg38 using \code{LOLA}
#' @param sigRegions A \code{GRanges} object of significant regions
#' @param regions A \code{GRanges} object of background regions
#' @param cores An integer of how many cores to use
#' @return A \code{tibble} of enrichment results
#' @import tidyverse
#' @import LOLA
#' @import simpleCache
#' @import magrittr
#' @import qvalue
#' @import hablar
#' @export chromHMM
chromHMM <- function(sigRegions = sigRegions,
                     regions = regions,
                     cores = cores){
  message("Performing ChromHMM enrichment testing")
  chromHMM <- LOLA::loadRegionDB(dbLocation = "/share/lasallelab/programs/LOLA/hg38",
                                 useCache = TRUE,
                                 limit = NULL,
                                 collections = "Roadmap_ChromHMM") %>%
    LOLA::runLOLA(userSets = sigRegions,
                  userUniverse = regions,
                  regionDB = .,
                  minOverlap = 1,
                  cores = cores,
                  redefineUserSets = FALSE) %T>%
    LOLA::writeCombinedEnrichment(combinedResults = .,
                                  outFolder = "ChromHMM",
                                  includeSplits = FALSE) %>%
    dplyr::as_tibble() %>%
    dplyr::select(oddsRatio, cellType, tissue, antibody) %>%
    dplyr::mutate(antibody = as.factor(antibody)) %>% 
    dplyr::mutate(antibody = dplyr::recode_factor(antibody,
                                                  "01_TssA" = "Active TSS",
                                                  "02_TssAFlnk" = "Flanking Active TSS",
                                                  "03_TxFlnk" = "Transcription at Gene 5' and 3'",
                                                  "04_Tx" = "Strong Transcription",
                                                  "05_TxWk" = "Weak Transcription",
                                                  "06_EnhG"= "Genic Enhancers",
                                                  "07_Enh" = "Enhancers",
                                                  "08_ZnfRpts" = "ZNF Genes & Repeats",
                                                  "09_Het" = "Heterochromatin",
                                                  "10_TssBiv" = "Bivalent/Poised TSS",
                                                  "11_BivFlnk" = "Flanking Bivalent TSS/Enhancer",
                                                  "12_EnhBiv" = "Bivalent Enhancer",
                                                  "13_ReprPC" = "Repressed PolyComb",
                                                  "14_ReprPCwk" = "Weak Repressed PolyComb",
                                                  "15_Quies" = "Quiescent/Low"
                                                  )
                  ) %>%
    dplyr::arrange(antibody) 
  
  # Fix Inf Odds Ratio
  max <- chromHMM %>%
  dplyr::summarize(max(hablar::s(oddsRatio))) %>% 
    dplyr::pull()
  
  chromHMM <- chromHMM %>%
    dplyr::mutate_if(is.numeric, function(x) ifelse(is.infinite(x), max, x)) %>% 
    tidyr::pivot_wider(names_from = antibody, values_from = oddsRatio) %>%
    dplyr::arrange(tissue) %>%
    return()
}

#' chromHMM_heatmap
#' @description Plot a heatmap of \code{LOLA} enrichment testing results of the ChromHMM 15-state model for hg38
#' @param chromHMM A \code{tibble} of enrichment results
#' @return Saves a heatmap
#' @import tidyverse
#' @import gplots
#' @import RColorBrewer
#' @import PerformanceAnalytics
#' @export chromHMM_heatmap
chromHMM_heatmap <- function(chromHMM = chromHMM){
  message("Plotting ChromHMM heatmap")
  
  # Make Row Labels
  labels <- chromHMM %>%
    dplyr::group_by(tissue) %>% 
    dplyr::tally()
  
  # Colors
  palette(PerformanceAnalytics::tol21rainbow)
  rowcolors <- as.list(palette(PerformanceAnalytics::tol21rainbow))
  rowcolors <- rowcolors[1:nrow(labels)]
  
  colorlist <- list()
  for (i in 1:length(labels$n)){
    colorlist[[i]] <- rep(rowcolors[i], labels$n[i])
  }
  colorlist <- unlist(colorlist)
  
  # Select matrix data
  data <- chromHMM %>%
    dplyr::select(-cellType, -tissue) %>%
    data.matrix()
  
  # Plot Heatmap
  pdf("ChromHMM/ChromHMM_heatmap.pdf",
      height = 8.5,
      width = 12)
  
  gplots::heatmap.2(data,
                    Rowv = F,
                    Colv = F,
                    dendrogram = "none",
                    col = (RColorBrewer::brewer.pal(9, name = "PuRd")),
                    margins = c(15,2),
                    trace = "none",
                    labRow = "" ,
                    labCol = colnames(data),
                    main = "Enriched Chromatin States for Differentially Methylated Regions",
                    RowSideColors = colorlist,
                    srtCol = 60,
                    keysize = 0.85,
                    key.par = list(cex = 0.5),
                    key.xlab = "Odds Ratio",
                    key.ylab = "Frequency",
                    key.title = ""
  )
  
  # Legend
  par(xpd = TRUE, mar = par()$mar + c(0,6,0,0))
  legend(x = -0.075,
         y= 0.9,
         legend = labels$tissue,
         col = unlist(rowcolors),
         lty = 1,
         lwd = 6,
         cex = 1,
         bty = "n")
  dev.off()
}

#' roadmap
#' @description Perfom enrichment testing against the Roadmap Epigenomics core histone modifications
#'  (5 marks, 127 epigenomes) for hg38 using \code{LOLA}
#' @param sigRegions A \code{GRanges} object of significant regions
#' @param regions A \code{GRanges} object of background regions
#' @param cores An integer of how many cores to use
#' @return A \code{tibble} of enrichment results
#' @import tidyverse
#' @import LOLA
#' @import simpleCache
#' @import magrittr
#' @import qvalue
#' @import hablar
#' @export roadmap
roadmap <- function(sigRegions = sigRegions,
                    regions = regions,
                    cores = cores){
  message("Performing Roadmap epigenomics enrichment testing")
  roadmap <- LOLA::loadRegionDB(dbLocation = "/share/lasallelab/programs/LOLA/hg38",
                                useCache = TRUE, 
                                limit = NULL,
                                collections = "roadmap_epigenomics") %>%
    LOLA::runLOLA(userSets = sigRegions,
                  userUniverse = regions,
                  regionDB = .,
                  minOverlap = 1,
                  cores = cores,
                  redefineUserSets = FALSE) %T>%
    LOLA::writeCombinedEnrichment(combinedResults = .,
                                  outFolder = "RoadmapEpigenomics",
                                  includeSplits = FALSE) %>%
    dplyr::as_tibble() %>%
    dplyr::left_join(readr::read_tsv("/share/lasallelab/programs/LOLA/hg38/roadmap_epigenomics/index.txt"),
                     by = "filename") %>%
    dplyr::select(oddsRatio, antibody.x, donor, anatomy) %>%
    dplyr::rename(antibody = antibody.x) %>% 
    dplyr::mutate(antibody = tidyr::replace_na(antibody, "DNase"))
  
  # Subset for Core Histone Mods from same 127 samples: H3K4me1, H3K4me3, H3K27me3, H3K36me3, H3K9me3
  core <- roadmap %>%
    dplyr::group_by(antibody) %>% 
    dplyr::tally() %>%
    dplyr::filter(n == 127) %>%
    dplyr::pull(antibody)
  
  roadmap <- roadmap %>% 
    dplyr::filter(antibody %in% core) %>%
    dplyr::arrange(anatomy)
    
  # Fix Inf Odds Ratio
  max <- roadmap %>%
    dplyr::summarize(max(hablar::s(oddsRatio))) %>% 
    dplyr::pull()
  
  roadmap <- roadmap %>%
    dplyr::mutate_if(is.numeric, function(x) ifelse(is.infinite(x), max, x)) %>% 
    tidyr::pivot_wider(names_from = antibody, values_from = oddsRatio) %>%
    dplyr::mutate(anatomy = stringr::str_replace(anatomy, "_", " ")) %>%
    dplyr::mutate(anatomy = stringr::str_to_title(anatomy)) %>%
    dplyr::mutate(anatomy = stringr::str_replace(anatomy, "Ipsc", "IPSC")) %>%
    dplyr::mutate(anatomy = stringr::str_replace(anatomy, "Esc", "ESC")) %>%
    dplyr::mutate(anatomy = stringr::str_replace(anatomy, "Gi", "GI")) %>%
    return()
}

#' roadmap_heatmap
#' @description Plot a heatmap of \code{LOLA} enrichment testing results of the Roadmap Epigenomocs
#'  core marks for hg38
#' @param chromHMM A \code{tibble} of enrichment results
#' @return Saves a heatmap
#' @import tidyverse
#' @import gplots
#' @import RColorBrewer
#' @import PerformanceAnalytics
#' @export roadmap_heatmap
roadmap_heatmap <- function(roadmap = roadmap){
  message("Plotting Roadmap Epigenomics heatmap")

  # Make Row Labels
  labels <- roadmap %>%
    dplyr::group_by(anatomy) %>% 
    dplyr::tally() 
  
  # Colors
  palette(PerformanceAnalytics::tol21rainbow)
  rowcolors <- as.list(palette(PerformanceAnalytics::tol21rainbow))
  rowcolors <- rep(rowcolors, 2)
  rowcolors <-rowcolors[1:nrow(labels)]
  
  colorlist <- list()
  for (i in 1:length(labels$n)){
    colorlist[[i]] <- rep(rowcolors[i], labels$n[i])
  }
  colorlist <- unlist(colorlist)
  
  # Select matrix data
  data <- roadmap %>%
    dplyr::select(-donor, -anatomy) %>%
    data.matrix()
  
  # Plot Heatmap
  pdf("RoadmapEpigenomics/Roadmap_heatmap.pdf",
      height = 8.5,
      width = 12)
  
  gplots::heatmap.2(data,
                    Rowv = F,
                    Colv = F,
                    dendrogram = "none",
                    col = (RColorBrewer::brewer.pal(9,name="PuRd")),
                    margins =c (15,2),
                    trace = "none",
                    labRow = "" ,
                    labCol = colnames(data),
                    main = "Histone Modifications at Differentially Methylated Regions",
                    RowSideColors = colorlist,
                    srtCol = 60,
                    keysize = 0.85,
                    key.par = list(cex = 0.5),
                    key.xlab = "Odds Ratio",
                    key.ylab = "Frequency",
                    key.title = ""
  )
  
  # Legend
  par(xpd = TRUE, mar = par()$mar + c(0,6,0,0))
  legend(x = -0.075,
         y= 0.9,
         legend = unlist(labels$anatomy),
         col = unlist(rowcolors),
         lty = 1,
         lwd = 6,
         cex = 1,
         bty = "n")
  dev.off()
}

#' dmrList
#' @description Create \code{GRangesList} object of all, hypermethylated, and hypomethylated DMRs from \code{dmrseq::dmrseq()}
#' @param sigRegions A \code{GRanges} object of DMRs from \code{dmrseq::dmrseq()}
#' @return A \code{GRangesList} of DMRs
#' @import tidyverse
#' @import GenomicRanges
#' @export dmrList
dmrList <- function(sigRegions = sigRegions){
  message("Making DMR list")
  
  sigRegions <- sigRegions %>%
    dplyr::as_tibble() %>%
    dplyr::mutate(direction = dplyr::case_when(stat > 0 ~ "Hypermethylated",
                                               stat < 0 ~ "Hypomethylated")
                  ) %>%
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE) 
  
    GenomicRanges::GRangesList("All DMRs" = sigRegions,
                               "Hypermethylated DMRs" = sigRegions %>%
                                 dplyr::as_tibble() %>%
                                 dplyr::filter(direction == "Hypermethylated") %>%
                                 GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE),
                               "Hypomethylated DMRs" = sigRegions %>%
                                 dplyr::as_tibble() %>%
                                 dplyr::filter(direction == "Hypomethylated") %>%
                                 GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)
                               ) %>%
    return()
}
