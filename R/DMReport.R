#' DMReport
#' @description Create an html report of significant regions from \code{dmrseq}
#' @param sigRegions \code{GRanges} object of signficant regions (DMRs or blocks) from \code{dmrseq} that 
#' were annotated by \code{DMRichR::annotateRegions}
#' @param regions \code{GRanges} object of background regions from \code{dmrseq}
#' @param bsseq Smoothed \code{bsseq} object
#' @param coverage Numeric of coverage samples were filtered for
#' @param name Character for html report name
#' @return Saves an html report of DMRs with genic annotations
#' @import gt
#' @import tidyverse
#' @importFrom glue glue
#' @export DMReport
DMReport <- function(sigRegions = sigRegions,
                     regions = regions,
                     bsseq = bs.filtered.bsseq,
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
    gt::as_raw_html(inline_css = TRUE) %>%
    write(glue::glue("{name}.html"))
  cat("Done", "\n")
}
