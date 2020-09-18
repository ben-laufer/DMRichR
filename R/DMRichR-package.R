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
  message("Loading DMRichR")
}