#' PCA
#' @description Performs and plots a PCA from individual smoothed methylation values
#' @param matrix Matrix of transposed individual methylation values
#' @param title Character string of title for plot and pdf
#' @return A \code{ggplot} object that can be viewed by calling it,
#'  saved with \code{ggplot2::ggsave()}, or further modified by adding \code{ggplot2} syntax.
#' @import ggbiplot
#' @importFrom dplyr case_when
#' @importFrom forcats fct_rev
#' @importFrom glue glue
#' @export PCA
PCA <- function(matrix = matrix,
                group = NA,
                title = title){
  print(glue::glue("Performing PCA..."))
  data.pca <- prcomp(matrix, center = TRUE, scale. = TRUE)
  #plot(data.pca, type = "l")
  #print(summary(data.pca))
  group <- factor(group, levels = unique(forcats::fct_rev(group)))
  ellipse <- dplyr::case_when(length(group) < 6 ~ F,
                              length(group) >= 6 ~ T)
  
  cat("Plotting PCA...")
  PCA <- ggbiplot::ggbiplot(data.pca,
                            obs.scale = 1,
                            var.scale = 1,
                            groups = group,
                            ellipse = ellipse,
                            circle = FALSE,
                            var.axes = FALSE,
                            choices = 1:2) +
    scale_color_discrete(name = '') +
    theme_bw(base_size = 20) +
    geom_point(aes(colour = group), size = 4) +
    ggtitle(title) + # Change title
    theme(legend.direction = 'vertical',
          #legend.position = c(0.125, 0.1), # Change legend position
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 18),
          panel.grid.major = element_blank(),
          panel.border = element_rect(color = "black", size = 1.25),
          axis.ticks = element_line(size = 1.25),
          legend.key = element_blank(),
          panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5)
          ) +
    guides(col = guide_legend(ncol = 1)) +
    theme()
  cat("Done", "\n")
  
  return(PCA)
}

#' windowsPCA
#' @description Performs and plots a PCA of 20kb windows from individual smoothed methylation values
#' @param goi A \code{BSgenome} object of the genome of interest (i.e. "BSgenome.Hsapiens.UCSC.hg38")
#' @param bsseq Smoothed \code{bsseq} object with a testCovariate in \code{pData}
#' @return A \code{ggplot} object that can be viewed by calling it,
#'  saved with \code{ggplot2::ggsave()}, or further modified by adding \code{ggplot2} syntax.
#' @import ggbiplot
#' @importFrom magrittr %>%
#' @importFrom dplyr select
#' @importFrom GenomeInfoDb seqlengths keepStandardChromosomes
#' @importFrom GenomicRanges tileGenome
#' @importFrom bsseq getMeth
#' @importClassesFrom bsseq BSseq 
#' @importMethodsFrom bsseq pData
#' @export windowsPCA
windowsPCA <- function(bsseq = bs.filtered.bsseq,
                       goi = goi,
                       group = NA){
  print(glue::glue("[DMRichR] Creating and plotting PCA of 20 kb windows from the {BSgenome::commonName(goi)} genome"))
  goi %>%
    GenomeInfoDb::seqlengths() %>%
    GenomicRanges::tileGenome(tilewidth = 2e4,
                              cut.last.tile.in.chrom = TRUE) %>%
    GenomeInfoDb::keepStandardChromosomes(pruning.mode = "coarse") %>%
    cbind(., data.frame(
      bsseq::getMeth(BSseq = bsseq,
                     regions = .,
                     type = "smooth",
                     what = "perRegion"),
      check.names = FALSE)
    ) %>%
    dplyr::select(-seqnames, -start, -end, -width, -strand) %>% 
    na.omit() %>%
    as.matrix() %>%
    t() %>% 
    DMRichR::PCA(group = group,
                 title = "Smoothed 20 Kb CpG Windows with CpG Islands") %>%
    return()
}

#' CGiPCA
#' @description Performs and plots a PCA of CpG island from individual smoothed methylation values
#' @param genome A character vector of the genome of interest (i.e. "hg38")
#' @param bsseq Smoothed \code{bsseq} object with a testCovariate in \code{pData}
#' @return A \code{ggplot} object that can be viewed by calling it,
#'  saved with \code{ggplot2::ggsave()}, or further modified by adding \code{ggplot2} syntax.
#' @import ggbiplot
#' @importFrom magrittr %>%
#' @importFrom dplyr select
#' @importFrom annotatr build_annotations
#' @importFrom GenomeInfoDb keepStandardChromosomes
#' @importFrom bsseq getMeth
#' @importClassesFrom bsseq BSseq 
#' @importMethodsFrom bsseq pData
#' @export CGiPCA
CGiPCA <- function(bsseq = bs.filtered.bsseq,
                   genome = genome,
                   group = NA){
  stopifnot(genome == "hg38" | genome == "mm10" | genome == "rn6")
  print(glue::glue("[DMRichR] Creating and plotting PCA of CpG islands from {genome}"))
  annotatr::build_annotations(genome = genome,
                              annotations = paste(genome,"_cpg_islands", sep = "")) %>% 
    GenomeInfoDb::keepStandardChromosomes(pruning.mode = "coarse") %>% 
    cbind(., data.frame(
      bsseq::getMeth(BSseq = bsseq,
                     regions = .,
                     type = "smooth",
                     what = "perRegion"),
      check.names = FALSE)
    ) %>% 
    dplyr::select(-seqnames, -start, -end, -width, -strand,
                  - id, -tx_id, -gene_id, -symbol, - type) %>% 
    na.omit() %>%
    as.matrix() %>%
    t() %>% 
    DMRichR::PCA(group = group,
                 title = "Smoothed CpG Island Windows") %>%
    return()
}

#' densityPlot
#' @description Creates a density plot of the mean of individual smoothed methylation values for every CpG
#' @param bsseq Smoothed \code{bsseq} object with a testCovariate in \code{pData}
#' @return A \code{ggplot} object that can be viewed by calling it,
#'  saved with \code{ggplot2::ggsave()}, or further modified by adding \code{ggplot2} syntax.
#' @importFrom magrittr %>% set_colnames
#' @importFrom dplyr as_tibble select transmute contains mutate
#' @importFrom forcats fct_rev
#' @import ggplot2
#' @importFrom tidyr gather
#' @importFrom bsseq getMeth
#' @importClassesFrom bsseq BSseq 
#' @importMethodsFrom bsseq pData
#' @export densityPlot

densityPlot <- function(bsseq = bs.filtered.bsseq,
                        group = NA){
  print(glue::glue("[DMRichR] Density plot of {nrow(bsseq)} CpGs"))
  bs.filtered.bsseq %>% 
    bsseq::getMeth(BSseq = .,
                   type = "smooth",
                   what = "perBase") %>%
    dplyr::as_tibble() %>% 
    na.omit() %>%
    magrittr::set_colnames(paste(group, seq_along(1:length(group)))) %>%
    dplyr::transmute(Group1 = dplyr::select(., dplyr::contains(levels(group)[1])) %>% rowMeans()*100,
                     Group2 = dplyr::select(., dplyr::contains(levels(group)[2])) %>% rowMeans()*100) %>%
    magrittr::set_colnames(c(levels(group)[1], levels(group)[2])) %>% 
    tidyr::gather(key = "variable",
                  value = "value") %>%
    dplyr::mutate(variable = factor(.$variable)) %>% 
    dplyr::mutate(variable = factor(.$variable, levels = unique(forcats::fct_rev(group)))) %>% 
    ggplot(aes(value, color = variable)) +
    geom_density(size = 1.2) +
    labs(x = "Percent Methylation", y = "Density", color = "Group") +
    theme_classic() +
    scale_x_continuous(expand = c(0.05,0.05), breaks = c(0,25,50,75,100)) +
    scale_y_continuous(expand = c(0.00,0.001)) +
    theme(axis.text = element_text(size = 16), axis.title = element_text(size = 16),
          strip.text = element_text(size = 16), legend.text = element_text(size = 14),
          legend.position = "bottom", legend.title = element_text(size = 14)) %>%
    return()
}
