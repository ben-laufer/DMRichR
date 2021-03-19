#' windows
#' @description Obtain windows of individual smoothed methylation values
#' @param bs.filtered.bsseq Smoothed \code{bsseq} object with a testCovariate in \code{pData}
#' @param size The number of bases in the window (default is 2e4, which is 20 Kb)
#' @param goi A \code{BSgenome} object of the genome of interest (i.e. "BSgenome.Hsapiens.UCSC.hg38")
#' @return A matrix of smoothed individual methylation values
#' @importFrom magrittr %>%
#' @importFrom GenomeInfoDb seqlengths keepStandardChromosomes
#' @importFrom GenomicRanges tileGenome
#' @importFrom bsseq getMeth
#' @importClassesFrom bsseq BSseq 
#' @importMethodsFrom bsseq pData
#' @export windows
#' 
windows <- function(bs.filtered.bsseq = bs.filtered.bsseq,
                    size = 2e4,
                    goi = goi){
  print(glue::glue("Obtaining {size/1000} Kb window individual smoothed methylation values from the {BSgenome::commonName(goi)} genome"))
  goi %>%
    GenomeInfoDb::seqlengths() %>%
    GenomicRanges::tileGenome(tilewidth = size,
                              cut.last.tile.in.chrom = TRUE) %>%
    GenomeInfoDb::keepStandardChromosomes(pruning.mode = "coarse") %>%
    bsseq::getMeth(BSseq = bs.filtered.bsseq,
                   regions = .,
                   type = "smooth",
                   what = "perRegion") %>% 
    na.omit() %>%
    return()
}

#' CGi
#' @description Obtain individual smoothed methylation values for CpG islands
#' @param bs.filtered.bsseq Smoothed \code{bsseq} object with a testCovariate in \code{pData}
#' @param genome A character vector of the genome of interest (i.e. "hg38")
#' @return A matrix of smoothed individual methylation values
#' @importFrom magrittr %>%
#' @importFrom GenomeInfoDb keepStandardChromosomes
#' @importFrom bsseq getMeth
#' @importClassesFrom bsseq BSseq 
#' @importMethodsFrom bsseq pData
#' @export CGi
#' 
CGi <- function(bs.filtered.bsseq = bs.filtered.bsseq,
                genome = genome){
  
  print(glue::glue("Obtaining individual smoothed methylation values of CpG islands from {genome}"))
  genome %>%
    DMRichR::getCpGs() %>% 
    bsseq::getMeth(BSseq = bs.filtered.bsseq,
                   regions = .,
                   type = "smooth",
                   what = "perRegion") %>% 
    na.omit() %>%
    return()
}

#' CpGs
#' @description Performs and plots a PCA of single CpG individual smoothed methylation values
#' @param bs.filtered.bsseq Smoothed \code{bsseq} object with a testCovariate in \code{pData}
#' @return A matrix of smoothed individual methylation values
#' @import ggbiplot
#' @importFrom magrittr %>% set_colnames
#' @importFrom bsseq getMeth
#' @importClassesFrom bsseq BSseq 
#' @importMethodsFrom bsseq pData
#' @export CpGs
#' 
CpGs <- function(bs.filtered.bsseq = bs.filtered.bsseq){
  print(glue::glue("Obtaining smoothed methylation values for all covered CpGs"))
  bs.filtered.bsseq %>% 
    bsseq::getMeth(BSseq = .,
                   type = "smooth",
                   what = "perBase") %>%
    na.omit() %>%
    return()
}

#' PCA
#' @description Performs and plots a PCA from individual smoothed methylation values
#' @param matrix A matrix of smoothed individual methylation values
#' @param group Ordered factor vector of sample groupings
#' @return A \code{ggplot} object that can be viewed by calling it,
#'  saved with \code{ggplot2::ggsave()}, or further modified by adding \code{ggplot2} syntax.
#' @import ggbiplot
#' @importFrom dplyr case_when
#' @importFrom forcats fct_rev
#' @importFrom glue glue
#' @importFrom Glimma glMDSPlot
#' @references \url{https://stackoverflow.com/questions/40315227/how-to-solve-prcomp-default-cannot-rescale-a-constant-zero-column-to-unit-var/40317343}
#' @export PCA
#' 
PCA <- function(matrix = matrix,
                group = NA){
  
  print(glue::glue("PCA of {length(matrix)} sites"))
  
  matrix %>%
    t() %>%
    .[ , which(apply(., 2, var) != 0)] %>% 
    prcomp(.,
           center = TRUE,
           scale. = TRUE) %>%
    ggbiplot::ggbiplot(.,
                       obs.scale = 1,
                       var.scale = 1,
                       groups = group,
                       ellipse = dplyr::case_when(length(group) < 6 ~ FALSE,
                                                  length(group) >= 6 ~ TRUE),
                       circle = FALSE,
                       var.axes = FALSE,
                       choices = 1:2) +
    scale_color_discrete(name = '') +
    theme_bw(base_size = 20) +
    geom_point(aes(colour = group), size = 8) +
    theme(legend.direction = 'vertical',
          #legend.position = c(0.125, 0.1), # Change legend position
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 18),
          panel.grid.major = element_blank(),
          panel.border = element_rect(color = "black", size = 1.25),
          axis.ticks = element_line(size = 1.25),
          legend.key = element_blank(),
          panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5)) +
    guides(col = guide_legend(ncol = 1)) +
    theme() %>%
    return()
}

#' densityPlot
#' @description Creates a density plot of the mean of individual smoothed methylation values
#' @param matrix A matrix of smoothed individual methylation values
#' @param group Ordered factor vector of sample groupings
#' @return A \code{ggplot} object that can be viewed by calling it,
#'  saved with \code{ggplot2::ggsave()}, or further modified by adding \code{ggplot2} syntax.
#' @importFrom magrittr %>% set_colnames
#' @importFrom dplyr as_tibble select transmute contains mutate
#' @importFrom forcats fct_rev
#' @import ggplot2
#' @importFrom tidyr gather
#' @export densityPlot
#' 
densityPlot <- function(matrix = matrix,
                        group = NA){
  
  print(glue::glue("Density plot of {length(matrix)} sites"))
  
  matrix  %>%
    dplyr::as_tibble() %>% 
    magrittr::set_colnames(paste(group, seq_along(1:length(group)))) %>%
    dplyr::transmute(Group1 = dplyr::select(., dplyr::contains(levels(group)[1])) %>% rowMeans()*100,
                     Group2 = dplyr::select(., dplyr::contains(levels(group)[2])) %>% rowMeans()*100) %>%
    magrittr::set_colnames(c(levels(group)[1], levels(group)[2])) %>% 
    tidyr::gather(key = "variable",
                  value = "value") %>%
    #dplyr::mutate(variable = factor(.$variable)) %>% 
    dplyr::mutate(variable = factor(.$variable, levels = levels(group))) %>% 
    ggplot(aes(value, color = variable)) +
    geom_density(size = 1.2) +
    labs(x = "Percent Methylation",
         y = "Density",
         color = "Group") +
    theme_classic() +
    scale_x_continuous(expand = c(0.05,0.05),
                       breaks = c(0,25,50,75,100)) +
    scale_y_continuous(expand = c(0.00,0.001)) +
    theme(axis.text = element_text(size = 16),
          axis.title = element_text(size = 16),
          strip.text = element_text(size = 16),
          legend.text = element_text(size = 14),
          legend.position = "bottom",
          legend.title = element_text(size = 14)) %>%
    return()
}
