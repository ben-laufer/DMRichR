#' annotateRegions
#' @description Annotate and tidy regions from \code{dmrseq::dmrseq()}
#' @param regions A \code{GRanges} object of DMRs, blocks, or background regions from \code{dmrseq::dmrseq()}
#' @param TxDb \code{TxDb} annotation package for genome of interest.
#' @param annoDb Character specifying \code{OrgDb} annotation package for species of interest.
#' @return A \code{tibble} of annotated regions
#' @import tidyverse
#' @import GenomicRanges
#' @import ChIPseeker
#' @export annotateRegions
annotateRegions <- function(regions = sigRegions,
                            TxDb = TxDb,
                            annoDb = annoDb){
regions %>% 
  dplyr::as_tibble() %>%
  dplyr::mutate(percentDifference = round(beta/pi *100)) %>%
  dplyr::mutate(fold = dplyr::case_when(stat > 0 ~ "Hypermethylated",
                                        stat < 0 ~ "Hypomethylated"
                                        )
                )%>%
  GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE) %>% 
  ChIPseeker::annotatePeak(TxDb = TxDb,
                           annoDb = annoDb,
                           overlap = "all"
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
