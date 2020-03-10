#' DMReport
#' @description Create an html report of a \code{ChIPseeker csAnno} peak object with genic annotations
#' @param tidySigRegionsAnno A \code{ChIPseeker csAnno} peak object of DMRs from \code{dmrseq::dmrseq()}
#'  that has been tidied with \code{tidyDMRs}
#' @param regions \code{GRanges} object of background regions
#' @param bsseq Smoothed \code{bsseq} object
#' @param coverage Numeric of coverage samples were filtered for
#' @return Saves an html report of DMRs with genic annotations
#' @import gt
#' @import tidyverse
#' @importFrom glue glue
#' @export DMReport
DMReport <- function(tidySigRegionsAnno = tidySigRegionsAnno,
                     regions = regions,
                     bsseq = bs.filtered.bsseq,
                     coverage = coverage){
  cat("\n","Preparing HTML report...")
  tidySigRegionsAnno %>%
    dplyr::select(-ENSEMBL, -betaCoefficient, -statistic) %>%
    dplyr::mutate(difference = difference/100) %>% 
    gt() %>%
    tab_header(
      title = glue::glue("{nrow(tidySigRegionsAnno)} Significant DMRs"),
      subtitle = glue::glue("{nrow(tidySigRegionsAnno)} Significant DMRs \\
                         {round(sum(tidySigRegionsAnno$statistic > 0) / nrow(tidySigRegionsAnno), digits = 2)*100}% hypermethylated, \\
                         {round(sum(tidySigRegionsAnno$statistic < 0) / nrow(tidySigRegionsAnno), digits = 2)*100}% hypomethylated \\
                         in {length(regions)} background regions \\
                         from {nrow(bs.filtered)} CpGs assayed at {coverage}x coverage")
    ) %>% 
    fmt_number(
      columns = vars("width", "CpGs"),
      decimals = 0
    ) %>% 
    fmt_scientific(
      columns = vars("p-value", "q-value"),
      decimals = 2
    ) %>%
    fmt_percent(
      columns = vars("difference"),
      drop_trailing_zeros = TRUE
    ) %>% 
    as_raw_html(inline_css = TRUE) %>%
    write("DMReport.html")
  cat("Done", "\n")
}
