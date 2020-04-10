#'  smoothPheatmap
#' @description Plot a heatmap of normalized individual smoothed methylation value z scores for selected regions (i.e. significant DMRs)
#' @param bsseq Smoothed \code{bsseq} object
#' @param regions \code{GRanges} object of regions to plot a heatmap for
#' @param testCovariate The factor tested for differences between groups
#' @param ... Additional arguments passed onto \code{pheatmap()}
#' @return Saves a pdf image of the heatmap in the DMR folder
#' @import bsseq
#' @import pheatmap
#' @import tidyverse
#' @importFrom glue glue
#' @references \url{https://davetang.org/muse/2018/05/15/making-a-heatmap-in-r-with-the-pheatmap-package/}
#' @export smoothPheatmap
smoothPheatmap <- function(regions = sigRegions,
                           bsseq = bs.filtered.bsseq,
                           testCovariate = testCovariate,
                           ...){
  cat("\n[DMRichR] DMR heatmap \t\t\t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
  
  getMeth(BSseq = bsseq,
          regions = regions,
          type = "smooth",
          what = "perRegion") %>% 
    as.matrix() %>%
    pheatmap::pheatmap(.,
                       scale = "row",
                       annotation_col =  pData(bsseq) %>%
                         as.data.frame() %>%
                         dplyr::select_if(~ nlevels(.) > 1),
                       color = RColorBrewer::brewer.pal(11,
                                                        name = "RdBu") %>%
                         rev(),
                       show_colnames = F,
                       #angle_col = 45,
                       border_color = "grey",
                       main = glue::glue("Z-Scores of {length(regions)} Differentially Methylated Regions"),
                       fontsize = 16,
                       filename = "./DMRs/heatmap.pdf",
                       width = 11,
                       height = 8.5,
                       annotation_colors = pData(bsseq) %>%
                         dplyr::as_tibble() %>%
                         dplyr::select(testCovariate) %>%
                         dplyr::distinct() %>%
                         dplyr::mutate_if(is.factor, as.character) %>%
                         dplyr::arrange(dplyr::desc(!!rlang::sym(testCovariate))) %>% 
                         t() %>%
                         rbind(DMRichR::gg_color_hue(2)) %>%
                         dplyr::as_tibble() %>%  
                         magrittr::set_colnames(dplyr::slice(., 1)) %>%
                         dplyr::slice(2) %>%
                         as.list() %>%
                         unlist() %>% 
                         list(testCovariate = .) %>%
                         setNames(testCovariate),
                       ...
                       ) %>%
    return()
}


# factors = bs.filtered.bsseq %>%
#   pData() %>%
#   as.data.frame() %>%
#   dplyr::rename(testCovariate = !!testCovariate) %>% 
#   dplyr::mutate(testCovariate = factor(.$testCovariate, levels = unique(forcats::fct_rev(.$testCovariate)))) %>%
#   dplyr::rename(!!testCovariate := testCovariate) %>%
#   dplyr::select_if(~ nlevels(.) > 1)

