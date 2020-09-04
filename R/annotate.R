#' annotateRegions
#' @description Annotate and tidy regions from \code{dmrseq::dmrseq()}
#' @param regions A \code{GRanges} object of DMRs, blocks, or background regions from \code{dmrseq::dmrseq()}
#' @param TxDb \code{TxDb} annotation package for genome of interest
#' @param annoDb Character specifying \code{OrgDb} annotation package for species of interest
#' @return A \code{tibble} of annotated regions
#' @importFrom dplyr rename as_tibble case_when mutate select
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom ChIPseeker annotatePeak
#' @importFrom magrittr %>%
#' @export annotateRegions
#' 
annotateRegions <- function(regions = sigRegions,
                            TxDb = TxDb,
                            annoDb = annoDb){
  regions %>% 
    dplyr::as_tibble() %>%
    dplyr::mutate(percentDifference = round(beta/pi *100)) %>%
    dplyr::mutate(direction = dplyr::case_when(stat > 0 ~ "Hypermethylated",
                                               stat < 0 ~ "Hypomethylated"
                                               )
                  )%>%
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE) %>% 
    ChIPseeker::annotatePeak(TxDb = TxDb,
                             annoDb = annoDb,
                             overlap = "all",
                             verbose = FALSE,
                             ) %>%
    dplyr::as_tibble() %>%
    dplyr::select("seqnames",
                  "start",
                  "end",
                  "width",
                  "L",
                  "beta",
                  "stat",
                  "pval",
                  "qval",
                  "percentDifference",
                  "annotation",
                  "distanceToTSS",
                  "ENSEMBL",
                  "SYMBOL",
                  "GENENAME"
                  ) %>%
    dplyr::rename(CpGs = L,
                  betaCoefficient = beta,
                  statistic = stat,
                  "p-value" = pval,
                  "q-value" = qval,
                  difference = percentDifference,
                  geneSymbol = SYMBOL,
                  gene = GENENAME
                  ) %>%
    return()
}

#' DMReport
#' @description Create an html report of significant regions from \code{dmrseq}
#' @param sigRegions \code{GRanges} object of signficant regions (DMRs or blocks) from \code{dmrseq} that 
#' were annotated by \code{DMRichR::annotateRegions}
#' @param regions \code{GRanges} object of background regions from \code{dmrseq}
#' @param bs.filtered Filtered \code{bsseq} object from \code{processBismark()}
#' @param coverage Numeric of coverage samples were filtered for
#' @param name Character for html report name
#' @return Saves an html report of DMRs with genic annotations
#' @importFrom gt gt tab_header fmt_number fmt_scientific fmt_percent as_raw_html
#' @importFrom dplyr select mutate 
#' @importFrom glue glue
#' @importFrom magrittr %>%
#' @importClassesFrom bsseq BSseq 
#' @export DMReport
#' 
DMReport <- function(sigRegions = sigRegions,
                     regions = regions,
                     bs.filtered = bs.filtered,
                     coverage = coverage,
                     name = "DMReport"){
  cat("\n","Preparing HTML report...")
  sigRegions %>%
    dplyr::select(-ENSEMBL, -betaCoefficient, -statistic) %>%
    dplyr::mutate(difference = difference/100) %>% 
    gt::gt() %>%
    gt::tab_header(
      title = glue::glue("{nrow(sigRegions)} Significant regions"),
      subtitle = glue::glue("{nrow(sigRegions)} Significant regions \\
                         {round(sum(sigRegions$statistic > 0) / nrow(sigRegions), digits = 2)*100}% hypermethylated, \\
                         {round(sum(sigRegions$statistic < 0) / nrow(sigRegions), digits = 2)*100}% hypomethylated \\
                         in {length(regions)} background regions \\
                         from {nrow(bs.filtered)} CpGs assayed at {coverage}x coverage")
    ) %>% 
    gt::fmt_number(
      columns = gt::vars("width", "CpGs"),
      decimals = 0
    ) %>% 
    gt::fmt_scientific(
      columns = vars("p-value", "q-value"),
      decimals = 2
    ) %>%
    gt::fmt_percent(
      columns = vars("difference"),
      drop_trailing_zeros = TRUE
    ) %>% 
    gt::as_raw_html(inline_css = FALSE) %>%
    write(glue::glue("{name}.html"))
  cat("Done", "\n")
}

#' annotateCpGs
#' @description Annotates DMRs from \code{dmrseq::dmrseq()} with CpG annotations
#'  using \code{annotatr} and returns a \code{ggplot2}
#' @param siRegions A \code{GRanges} object of significant DMRs returned by \code{dmrseq:dmrseq()}
#' @param regions A \code{GRanges} object of background regions returned by \code{dmrseq:dmrseq()}
#' @param genome A character vector specifying the genome of interest
#'  c("hg38", "hg19", "mm10", "mm9", "rn6", "rn5")
#' @param saveAnnotations A logical indicating whether to save bed files of annotations
#'  for external enrichment testing
#' @return A \code{ggplot} object of CpG annotations that can be viewed by calling it,
#'  saved with \code{ggplot2::ggsave()}, or further modified by adding \code{ggplot2} syntax.
#' @import ggplot2
#' @importFrom annotatr build_annotations annotate_regions plot_categorical
#' @importFrom GenomeInfoDb keepStandardChromosomes
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom dplyr as_tibble mutate case_when
#' @importFrom glue glue
#' @importFrom magrittr %>%
#' @export annotateCpGs
#' 
annotateCpGs <- function(sigRegions = sigRegions,
                         regions = regions,
                         genome = genome,
                         saveAnnotations = F){
  stopifnot(genome %in% c("hg38", "hg19", "mm10", "mm9", "rn6", "rn5"))
  cat("\n[DMRichR] Building CpG annotations \t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
  annotations <- annotatr::build_annotations(genome = genome, annotations = paste(genome,"_cpgs", sep="")) %>%
    GenomeInfoDb::keepStandardChromosomes(pruning.mode = "coarse")
  
  glue::glue("Annotating DMRs...")
  sigRegions <- sigRegions %>%
    dplyr::as_tibble() %>%
    dplyr::mutate(direction = dplyr::case_when(stat > 0 ~ "Hypermethylated",
                                               stat < 0 ~ "Hypomethylated")
    ) %>%
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE) 
  
  dm_annotated_CpG <- annotatr::annotate_regions(
    regions = sigRegions,
    annotations = annotations,
    ignore.strand = TRUE,
    quiet = FALSE)
  
  glue::glue("Annotating background regions...")
  regions <- regions %>%
    dplyr::as_tibble() %>%
    dplyr::mutate(direction = dplyr::case_when(stat > 0 ~ "Hypermethylated",
                                               stat < 0 ~ "Hypomethylated")
    ) %>%
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE) 
  
  background_annotated_CpG <- annotatr::annotate_regions(
    regions = regions,
    annotations = annotations,
    ignore.strand = TRUE,
    quiet = FALSE)
  
  if(saveAnnotations == T){
    glue::glue("Saving files for GAT...")
    if(dir.exists("Extra") == F){dir.create("Extra")}
    CpGs <- as.data.frame(annotations)
    CpGs <- CpGs[!grepl("_", CpGs$seqnames), ]
    table(CpGs$seqnames)
    DMRichR::df2bed(CpGs[, c(1:3,10)], paste("Extra/GAT/", genome, "CpG.bed", sep = ""))
  }

  glue::glue("Preparing CpG annotation plot...")
  CpG_bar <- annotatr::plot_categorical(
    annotated_regions = dm_annotated_CpG,
    annotated_random = background_annotated_CpG,
    x = 'direction',
    fill = 'annot.type',
    x_order = c('Hypermethylated','Hypomethylated'),
    fill_order = c(
      paste(genome,"_cpg_islands", sep = ""),
      paste(genome,"_cpg_shores", sep = ""),
      paste(genome,"_cpg_shelves", sep = ""),
      paste(genome,"_cpg_inter", sep = "")
      ),
    position = 'fill',
    plot_title = '',
    legend_title = 'Annotations',
    x_label = '',
    y_label = 'Proportion') +
    scale_x_discrete(labels = c("All", "Hypermethylated", "Hypomethylated", "Background")) +
    scale_y_continuous(expand = c(0,0)) +
    theme_classic() +
    theme(axis.text = element_text(size = 25),
          axis.title = element_text(size = 25),
          strip.text = element_text(size = 25),
          #legend.position = "none",
          axis.text.x = element_text(angle = 45,
                                     hjust = 1)) %>%
    return()
}

#' annotateGenic
#' @description Annotates DMRs from \code{dmrseq::dmrseq()} with genic annotations
#'  using \code{annotatr} and returns a \code{ggplot2}
#' @param siRegions A \code{GRanges} object of signficant DMRs returned by \code{dmrseq:dmrseq()}
#' @param regions A \code{GRanges} object of background regions returned by \code{dmrseq:dmrseq()}
#' @param genome A character vector specifying the genome of interest ("hg38" or "mm10")
#' @param saveAnnotations A logical indicating whether to save bed files of annoation database for external enrichment testing
#' @return A \code{ggplot} object of genic annotations that can be viewed by calling it,
#'  saved with \code{ggplot2::ggsave()}, or further modified by adding \code{ggplot2} syntax.
#' @import ggplot2
#' @importFrom annotatr build_annotations annotate_regions plot_categorical
#' @importFrom GenomeInfoDb keepStandardChromosomes
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom dplyr as_tibble mutate case_when
#' @importFrom glue glue
#' @importFrom magrittr %>%
#' @export annotateGenic
#' 
annotateGenic <- function(sigRegions = sigRegions,
                          regions = regions,
                          genome = genome,
                          saveAnnotations = F){
  stopifnot(genome %in% c("hg19", "hg38", "mm9", "mm10", "rn5", "rn6", "dm6"))
  cat("\n[DMRichR] Building gene region annotations \t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
  annotations <- annotatr::build_annotations(genome = genome, annotations = c(paste(genome,"_basicgenes", sep = ""),
                                                                              paste(genome,"_genes_intergenic", sep = ""),
                                                                              paste(genome,"_genes_intronexonboundaries", sep = ""),
                                                                              if(genome == "hg38" | genome == "mm10"){paste(genome,"_enhancers_fantom", sep = "")})) %>%
    GenomeInfoDb::keepStandardChromosomes(., pruning.mode = "coarse")
  
  if(saveAnnotations == T){
    glue::glue("Saving files for GAT...")
    if(dir.exists("Extra") == F){dir.create("Extra")}
    annoFile <- as.data.frame(annotations)
    annoFile <- annoFile[!grepl("_", annoFile$seqnames) ,]
    table(annoFile$seqnames)
    annoFile <- annoFile[, c(1:3,10)]
    
    if(genome == "hg38" | genome == "mm10"){DMRichR::gr2bed(annoFile[annoFile$type == paste(genome,"_enhancers_fantom", sep = ""), ], "Extra/GAT/enhancers.bed")}
    DMRichR::gr2bed(annoFile[annoFile$type == paste(genome,"_genes_promoters", sep = ""), ], "Extra/GAT/promoters.bed")
    DMRichR::gr2bed(annoFile[annoFile$type == paste(genome,"_genes_introns", sep = ""), ], "Extra/GAT/introns.bed")
    DMRichR::gr2bed(annoFile[annoFile$type == paste(genome,"_genes_intronexonboundaries", sep = ""), ], "Extra/GAT/boundaries.bed")
    DMRichR::gr2bed(annoFile[annoFile$type == paste(genome,"_genes_intergenic", sep = ""), ], "Extra/GAT/intergenic.bed")
    DMRichR::gr2bed(annoFile[annoFile$type == paste(genome,"_genes_exons", sep = ""), ], "Extra/GAT/exons.bed")
    DMRichR::gr2bed(annoFile[annoFile$type == paste(genome,"_genes_5UTRs", sep = ""), ], "Extra/GAT/fiveUTRs.bed")
    DMRichR::gr2bed(annoFile[annoFile$type == paste(genome,"_genes_3UTRs", sep = ""), ], "Extra/GAT/threeUTRs.bed")
    DMRichR::gr2bed(annoFile[annoFile$type == paste(genome,"_genes_1to5kb", sep = ""), ], "Extra/GAT/onetofivekb.bed")
  }
  
  glue::glue("Annotating DMRs...")
  sigRegions <- sigRegions %>%
    dplyr::as_tibble() %>%
    dplyr::mutate(direction = dplyr::case_when(stat > 0 ~ "Hypermethylated",
                                               stat < 0 ~ "Hypomethylated")
    ) %>%
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE) 
  
  dm_annotated <- annotatr::annotate_regions(
    regions = sigRegions,
    annotations = annotations,
    ignore.strand = TRUE,
    quiet = FALSE)
  
  glue::glue("Annotating background regions...")
  regions <- regions %>%
    dplyr::as_tibble() %>%
    dplyr::mutate(direction = dplyr::case_when(stat > 0 ~ "Hypermethylated",
                                               stat < 0 ~ "Hypomethylated")
    ) %>%
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE) 
  
  background_annotated <- annotatr::annotate_regions(
    regions = regions,
    annotations = annotations,
    ignore.strand = TRUE,
    quiet = FALSE)
  
  glue::glue("Preparing CpG annotation plot...")
  gene_bar <- annotatr::plot_categorical(
    annotated_regions = dm_annotated,
    annotated_random = background_annotated,
    x = 'direction',
    fill = 'annot.type',
    x_order = c('Hypermethylated','Hypomethylated'),
    fill_order =  c(
      if(genome == "hg38" | genome == "mm10"){paste(genome, "_enhancers_fantom", sep = "")},
      paste(genome,"_genes_1to5kb", sep = ""),
      paste(genome,"_genes_promoters", sep = ""),
      paste(genome,"_genes_5UTRs", sep = ""),
      paste(genome,"_genes_exons", sep = ""),
      paste(genome,"_genes_intronexonboundaries", sep = ""),
      paste(genome,"_genes_introns", sep = ""),
      paste(genome,"_genes_3UTRs", sep = ""),
      paste(genome,"_genes_intergenic", sep = "")
      ),
    position = 'fill',
    plot_title = '',
    legend_title = 'Annotations',
    x_label = '',
    y_label = 'Proportion') +
    scale_x_discrete(labels = c("All", "Hypermethylated", "Hypomethylated", "Background")) +
    scale_y_continuous(expand = c(0,0)) +
    theme_classic() +
    theme(axis.text = element_text(size = 25),
          axis.title = element_text(size = 25),
          strip.text = element_text(size = 25),
          #legend.position = "none",
          axis.text.x = element_text(angle = 45,
                                     hjust = 1)) %>%
    return()
}
