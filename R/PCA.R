#' PCA
#' @description Performs and plots a PCA from individual smoothed methylation values
#' @param matrix Matrix of transposed individual methylation values
#' @param title Character string of title for plot and pdf
#' @return A \code{ggplot} object that can be viewed by calling it,
#'  saved with \code{ggplot2::ggsave()}, or further modified by adding \code{ggplot2} syntax.
#' @import ggbiplot
#' @export PCA
PCA <- function(matrix = matrix,
                group = bs.filtered.bsseq %>% pData() %>% dplyr::as_tibble() %>% dplyr::pull(!!testCovariate),
                title = title){
  print(glue::glue("Performing PCA..."))
  data.pca <- prcomp(matrix, center = TRUE, scale. = TRUE)
  plot(data.pca, type = "l")
  print(summary(data.pca))

  cat("Plotting PCA...")
  PCA <- ggbiplot::ggbiplot(data.pca,
                            obs.scale = 1,
                            var.scale = 1,
                            groups = group,
                            ellipse = TRUE,
                            circle = FALSE,
                            var.axes = FALSE,
                            choices = 1:2) +
    scale_color_discrete(name = '') +
    theme_bw(base_size = 25) +
    geom_point(aes(colour = group), size = 4) +
    theme(legend.direction = 'vertical',
          #legend.position = c(0.125, 0.1), # Change legend position
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 18),
          panel.grid.major = element_blank(),
          panel.border = element_rect(color = "black", size = 1.25),
          axis.ticks = element_line(size = 1.25),
          legend.key = element_blank(),
          panel.grid.minor = element_blank()) +
    guides(col=guide_legend(ncol=2)) +
    ggtitle(title) + # Change title
    theme(plot.title = element_text(hjust = 0.5))
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
#' @export windowsPCA
windowsPCA <- function(bsseq = bs.filtered.bsseq,
                       goi = goi){
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
    DMRichR::PCA(group = bsseq %>%
                   pData() %>%
                   dplyr::as_tibble() %>%
                   dplyr::pull(!!testCovariate),
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
#' @export CGiPCA
CGiPCA <- function(bsseq = bs.filtered.bsseq,
                   genome = genome){
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
    DMRichR::PCA(group = bsseq %>%
                   pData() %>%
                   dplyr::as_tibble() %>%
                   dplyr::pull(!!testCovariate),
                 title = "Smoothed CpG Island Windows") %>%
    return()
}
