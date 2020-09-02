#' DMRichGenic
#' @description Test DMRs for overlaps with gene region annotations for all genomes.
#' @param sigRegions \code{GRanges} object of DMRs.
#' @param regions \code{GRanges} object of background regions. 
#' @param annoDb Character specifying OrgDb annotation package for species of interest.
#' @param TxDb TxDb annotation package for genome of interest.
#' @return A tibble with the enrichment results.
#' @importFrom dplyr filter mutate case_when select as_tibble
#' @importFrom magrittr %>%
#' @importFrom glue glue glue_collapse
#' @importFrom stringr str_detect
#' @importFrom data.table rbindlist
#' @importFrom GenomeInfoDb genome
#' @importFrom forcats as_factor
#' @export DMRichGenic

DMRichGenic <- function(sigRegions = sigRegions,
                        regions = regions,
                        TxDb = TxDb,
                        annoDb = annoDb){
  
  genome <- TxDb %>%
    GenomeInfoDb::genome() %>%
    unique()
  
  print(glue::glue("{genome} annotations will be used for sigRegions and regions"))
  sigRegionsAnnotated <- sigRegions %>%
    DMRichR::annotateRegions(TxDb,
                             annoDb)
  
  regionsAnnotated <- regions %>%
    DMRichR::annotateRegions(TxDb,
                             annoDb)
  
  annotations <- c("Promoter", "5' UTR", "Exon", "Intron", "3' UTR", "Downstream", "Intergenic")
  print(glue::glue("Performing gene region enrichment testing for {genome} with the following annotations:"))
  print(glue::glue_collapse({annotations}, sep = ", "))
  
  lapply(annotations, function(term){
    
    print(glue::glue("Now FISHERing for {term} annotation"))
    
    sigRegionsOverlap <- sigRegionsAnnotated %>%
      dplyr::filter(stringr::str_detect(annotation, term))
                    
    regionsOverlap <- regionsAnnotated %>%
      dplyr::filter(stringr::str_detect(annotation, term))
    
    genicMatrix <- matrix(c(nrow(sigRegionsOverlap), (length(sigRegions) - nrow(sigRegionsOverlap)),
                              nrow(regionsOverlap), (length(regions) - nrow(regionsOverlap))),
                            nrow = 2)
    
    results <- data.frame("Annotation" =  term,
                          "OR" = fisher.test(genicMatrix)[["estimate"]][["odds ratio"]],
                          "CIlower" = fisher.test(genicMatrix)[["conf.int"]][1],
                          "CIupper" = fisher.test(genicMatrix)[["conf.int"]][2],
                          "p.value" = fisher.test(genicMatrix)$p.value)
    
  }) %>%
    data.table::rbindlist() %>%
    dplyr::mutate(fdr = p.adjust(p.value, method = 'fdr')) %>%
    dplyr::mutate(Annotation = forcats::as_factor(Annotation)) %>% 
    dplyr::as_tibble() %>% 
    return()
}

#' DMRichCpG
#' @description Test DMRs for overlaps with CpG annotations for human, mouse, or rat.
#' @param sigRegions \code{GRanges} object of DMRs.
#' @param regions \code{GRanges} object of background regions. 
#' @param genome A character vector specifying the genome of interest
#'  c("hg38", "hg19", "mm10", "mm9", "rn6", "rn5")
#' @return A tibble with the enrichment results.
#' @importFrom dplyr filter mutate case_when select recode_factor as_tibble
#' @importFrom magrittr %>%
#' @importFrom glue glue glue_collapse
#' @importFrom stringr str_detect
#' @importFrom data.table rbindlist
#' @importFrom annotatr build_annotations
#' @importFrom GenomeInfoDb keepStandardChromosomes genome
#' @importFrom plyranges as_granges mutate count_overlaps
#' @export DMRichCpG

DMRichCpG <- function(sigRegions = sigRegions,
                      regions = regions,
                      genome = genome){
  
  stopifnot(genome %in% c("hg38", "hg19", "mm10", "mm9", "rn6", "rn5"))
    
    print(glue::glue("Performing CpG annotation enrichment testing for {genome}"))
    CGannotations <- annotatr::build_annotations(genome = genome,
                                                 annotations = paste(genome,"_cpgs", sep = "")
    ) %>%
      GenomeInfoDb::keepStandardChromosomes(pruning.mode = "coarse") %>%
      plyranges::as_granges()
    
    CGterms <- c("cpg_islands", "cpg_shores", "cpg_shelves", "cpg_inter")
    
    lapply(CGterms, function(term){
      
      print(glue::glue("Now FISHERing for {term} annotation"))
      
      CGannotationsFiltered <- CGannotations %>%
        plyranges::filter(stringr::str_detect(type, term))
      
      sigRegionsOverlap <- sigRegions %>%
        plyranges::as_granges() %>% 
        plyranges::mutate(n_overlaps = plyranges::count_overlaps(.,CGannotationsFiltered)) %>%
        plyranges::filter(n_overlaps > 0)
      
      regionsOverlap <- regions %>%
        plyranges::as_granges() %>% 
        plyranges::mutate(n_overlaps = plyranges::count_overlaps(.,CGannotationsFiltered)) %>%
        plyranges::filter(n_overlaps > 0)
      
      CGmatrix <- matrix(c(length(sigRegionsOverlap), (length(sigRegions) - length(sigRegionsOverlap)),
                           length(regionsOverlap), (length(regions) - length(regionsOverlap))),
                         nrow = 2)
      
      results <- data.frame("Annotation" =  term,
                            "OR" = fisher.test(CGmatrix)[["estimate"]][["odds ratio"]],
                            "CIlower" = fisher.test(CGmatrix)[["conf.int"]][1],
                            "CIupper" = fisher.test(CGmatrix)[["conf.int"]][2],
                            "p.value" = fisher.test(CGmatrix)$p.value)
      
    }) %>%
      data.table::rbindlist() %>%
      dplyr::mutate(fdr = p.adjust(p.value, method = 'fdr')) %>%
      dplyr::mutate(Annotation = dplyr::recode_factor(Annotation,
                                                      cpg_islands = "CpG Islands",
                                                      cpg_shores = "CpG Shores",
                                                      cpg_shelves = "CpG Shelves",
                                                      cpg_inter = "Open Sea")
                    ) %>% 
      dplyr::as_tibble() %>% 
      return()
}

#' DMRichGenicPlot
#' @description Plot DMR gene region enrichment testing results from \code{DMRichR::DMRichGenic}.
#' @param genicOverlaps \code{tibble} from \code{DMRichR::DMRichGenic}.
#' @return A \code{ggplot} object of enrichment results that can be viewed by calling it, 
#' saved with \code{ggplot2::ggsave()}, or further modified by adding \code{ggplot2} syntax.
#' @import ggplot2
#' @import wesanderson
#' @importFrom dplyr filter mutate case_when select
#' @importFrom magrittr %>%
#' @importFrom glue glue
#' @export DMRichGenicPlot

DMRichGenicPlot <- function(genicOverlaps = genicOverlaps){
  
  print(glue::glue("Plotting genic enrichment results"))
  
  data <- genicOverlaps %>%
    dplyr::mutate(OR = dplyr::case_when(OR < 1 ~ -1/OR,
                                        OR >= 1 ~ OR)
    ) %>%
    dplyr::mutate(signif = dplyr::case_when(fdr <= 0.05 ~ 1,
                                            fdr> 0.05 ~ 0)
    ) %>%
    dplyr::select(Annotation,
                  OR,
                  fdr,
                  signif)
  
    ggplot(data = data,
          aes(x = Annotation,
              y = OR,
              fill = Annotation)
          ) +
      geom_bar(stat = "identity", 
               color = "Black") +
      coord_flip() +
      labs(y = "Fold Enrichment",
           x = element_blank()
           ) +
      theme_classic() + 
      theme(axis.text = element_text(size = 16),
            axis.title = element_text(size = 16),
            strip.text = element_text(size = 16),
            legend.text = element_text(size = 14),
            legend.title = element_text(size = 14)
            ) +
      scale_y_continuous(expand = c(0.1, 0.1)) + 
      scale_x_discrete(limits = data$Annotation %>%
                         levels() %>%
                         rev()) + 
      scale_fill_manual(values =  data$Annotation %>%
                          nlevels() %>% 
                          wesanderson::wes_palette("Zissou1", n = ., type = "continuous") %>%
                          rev(),
                        breaks = data$Annotation %>%
                          levels(),
                        name = "Annotation") +
      geom_hline(yintercept = 0) +
      geom_text(data = data[(data$signif == 1 & data$OR > 0), ],
                label = "*",
                size = 8,
                show.legend = FALSE,
                nudge_y = 0.5,
                nudge_x = -0.09) +
      geom_text(data = data[(data$signif == 1 & data$OR < 0), ],
                label = "*",
                size = 8,
                show.legend = FALSE,
                nudge_y = -0.5,
                nudge_x = -0.09)
  
}

#' DMRichCpGPlot
#' @description Plot DMR gene region enrichment testing results from \code{DMRichR::CpG}.
#' @param CGoverlaps \code{tibble} from \code{DMRichR::DMRichCpG}.
#' @return A \code{ggplot} object of enrichment results that can be viewed by calling it, 
#' saved with \code{ggplot2::ggsave()}, or further modified by adding \code{ggplot2} syntax.
#' @import ggplot2
#' @importFrom dplyr filter mutate case_when select
#' @importFrom magrittr %>%
#' @importFrom glue glue
#' @export DMRichCpGPlot

DMRichCpGPlot <- function(CGoverlaps = CGoverlaps){

  print(glue::glue("Plotting CpG enrichment results"))
  
  data <- CGoverlaps %>%
    dplyr::mutate(OR = dplyr::case_when(OR < 1 ~ -1/OR,
                                        OR >= 1 ~ OR)
    ) %>%
    dplyr::mutate(signif = dplyr::case_when(fdr <= 0.05 ~ 1,
                                            fdr> 0.05 ~ 0)
    ) %>%
    dplyr::select(Annotation,
                  OR,
                  fdr,
                  signif)
  
  ggplot(data = data,
         aes(x = Annotation,
             y = OR,
             fill = Annotation)
         ) +
    geom_bar(stat = "identity", 
             color = "Black") +
    coord_flip() +
    labs(y = "Fold Enrichment",
         x = element_blank()
         ) +
    theme_classic() + 
    theme(axis.text = element_text(size = 16),
          axis.title = element_text(size = 16),
          strip.text = element_text(size = 16),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 14)
          ) +
    scale_y_continuous(expand = c(0.1, 0.1)) + 
    scale_x_discrete(limits = data$Annotation %>%
                       levels() %>%
                       rev()) + 
    scale_fill_manual(values = c("forestgreen", "goldenrod2", "dodgerblue", "blue3"),
                      breaks = data$Annotation %>%
                        levels(),
                      name = "Annotation") +
    geom_hline(yintercept = 0) +
    geom_text(data = data[(data$signif == 1 & data$OR > 0), ],
              label = "*",
              size = 8,
              show.legend = FALSE,
              nudge_y = 0.5,
              nudge_x = -0.09) +
    geom_text(data = data[(data$signif == 1 & data$OR < 0), ],
              label = "*",
              size = 8,
              show.legend = FALSE,
              nudge_y = -0.5,
              nudge_x = -0.09)
  
}
