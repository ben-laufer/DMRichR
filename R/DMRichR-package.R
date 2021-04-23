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
  packageStartupMessage(paste0("DMRichR v", utils::packageVersion("DMRichR")," has loaded. ",
  "If you use DMRichR in published research please cite ",
  "Laufer et al. 2020, Korthauer et al. 2018, and Hansen et al. 2012."))
}