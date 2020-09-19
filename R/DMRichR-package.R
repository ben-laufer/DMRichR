#' DMRichR
#'
#' Enrich your DMR analysis
#' @name DMRichR
#' @docType package
NULL

#' onAttach
#' @title onAttach
#' @param libname Library name
#' @param pkgname Package name
#
.onAttach <- function(libname = find.package("DMRichR"), pkgname = "DMRichR") {
  message("DMRichR v1.5 has loaded.")
  message("Please cite Laufer et al. 2020 if you use this in published research.")
}