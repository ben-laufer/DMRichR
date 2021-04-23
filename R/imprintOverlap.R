#' imprintOverlap
#' @title Imprinted gene enrichment testing
#' @description Test DMRs for overlaps with human imprinted genes from
#'  \url{https://www.geneimprint.com/site/genes-by-species.Homo+sapiens.imprinted-All}
#' @param sigRegions \code{GRanges} object of DMRs.
#' @param regions \code{GRanges} object of background regions. 
#' @param annoDb Character specifying OrgDb annotation package for species of interest.
#' @param TxDb TxDb annotation package for genome of interest.
#' @return A printed statement with the enrichment results.
#' @importFrom dplyr filter mutate case_when select distinct pull
#' @importFrom magrittr %>%
#' @importFrom glue glue glue_collapse
#' @importFrom stringr str_detect
#' @importFrom plyranges as_granges
#' @importFrom stats fisher.test
#' @export imprintOverlap
#' 
imprintOverlap <- function(sigRegions = sigRegions,
                           regions = regions,
                           TxDb = TxDb,
                           annoDb = annoDb){
  
  imprint <- c("(?i)^DIRAS3$",
               "(?i)^RNU5D-1$",
               "(?i)^TP73$",
               "(?i)^LRRTM1$",
               "(?i)^GPR1$",
               "(?i)^ZDBF2$",
               "(?i)^NAP1L5$",
               "(?i)^ERAP2$",
               "(?i)^VTRNA2-1$",
               "(?i)^TNDM$",
               "(?i)^ADTRP$",
               "(?i)^FAM50B$",
               "(?i)^LIN28B$",
               "(?i)^AIM1$",
               "(?i)^PLAGL1$",
               "(?i)^HYMAI$",
               "(?i)^SLC22A2$",
               "(?i)^SLC22A3$",
               "(?i)^GRB10$",
               "(?i)^DDC$",
               "(?i)^PEG10$",
               "(?i)^MAGI2$",
               "(?i)^SGCE$",
               "(?i)^PPP1R9A$",
               "(?i)^TFPI2$",
               "(?i)^DLX5$",
               "(?i)^MEST$",
               "(?i)^CPA4$",
               "(?i)^COPG2IT1$",
               "(?i)^MESTIT1$",
               "(?i)^KLF14$",
               "(?i)^DLGAP2$",
               "(?i)^KCNK9$",
               "(?i)^ZFAT-AS1$",
               "(?i)^PEG13$",
               "(?i)^ZFAT$",
               "(?i)^GLIS3$",
               "(?i)^INPP5F$",
               "(?i)^WT1-AS$",
               "(?i)^WT1$",
               "(?i)^KCNQ1OT1$",
               "(?i)^OSBPL5$",
               "(?i)^KCNQ1DN$",
               "(?i)^KCNQ1$",
               "(?i)^H19$",
               "(?i)^CDKN1C$",
               "(?i)^PHLDA2$",
               "(?i)^IGF2AS$",
               "(?i)^SLC22A18$",
               "(?i)^INS$",
               "(?i)^IGF2$",
               "(?i)^ANO1$",
               "(?i)^ZC3H12C$",
               "(?i)^NTM$",
               "(?i)^RBP5$",
               "(?i)^LRP1$",
               "(?i)^ATP5F1EP2$",
               "(?i)^RB1$",
               "(?i)^SMOC1$",
               "(?i)^DLK1$",
               "(?i)^MEG3$",
               "(?i)^RTL1$",
               "(?i)^SNORD114",
               "(?i)^MEG8$",
               "(?i)^SNORD113",
               "(?i)^MAGEL2$",
               "(?i)^MKRN3$",
               "(?i)^UBE3A$",
               "(?i)^NPAP1$",
               "(?i)^SNORD109A$",
               "(?i)^SNORD108$",
               "(?i)^SNORD107$",
               "(?i)^SNORD109B$",
               "(?i)^ATP10A$",
               "(?i)^SNRPN$",
               "(?i)^SNORD116",
               "(?i)^SNORD115",
               "(?i)^PWCR1$",
               "(?i)^NDN$",
               "(?i)^SNURF$",
               "(?i)^SNORD64$",
               "(?i)^IRAIN$",
               "(?i)^ZNF597$",
               "(?i)^NAA60$",
               "(?i)^TP53$",
               "(?i)^TCEB3C$",
               "(?i)^DNMT1$",
               "(?i)^ZIM2$",
               "(?i)^PEG3$",
               "(?i)^MIMT1$",
               "(?i)^NLRP2$",
               "(?i)^MIR371A$",
               "(?i)^PSIMCT-1$",
               "(?i)^BLCAP$",
               "(?i)^NNAT$",
               "(?i)^MCTS2$",
               "(?i)^GDAP1L1$",
               "(?i)^SGK2$",
               "(?i)^GNAS$",
               "(?i)^L3MBTL$",
               "(?i)^SANG$",
               "(?i)^GNASAS$",
               "(?i)^MIR298$",
               "(?i)^MIR296$",
               "(?i)^DGCR6$",
               "(?i)^DGCR6L$")
  
  sigRegions <- sigRegions %>%
    plyranges::as_granges() %>% 
    DMRichR::annotateRegions(TxDb,
                             annoDb)
  
  regions <- regions %>%
    DMRichR::annotateRegions(TxDb,
                             annoDb)
  
  # Overlap pipe
  imprintOverlaps <- . %>%
    dplyr::filter(stringr::str_detect(geneSymbol, paste(imprint, collapse = "|")))
  
  cleanOverlaps <- . %>% 
    dplyr::mutate(geneSymbol = dplyr::case_when(stringr::str_detect(geneSymbol,"SNORD113") ~ "SNORD113",
                                                stringr::str_detect(geneSymbol,"SNORD114") ~ "SNORD114",
                                                stringr::str_detect(geneSymbol,"SNORD115") ~ "SNORD115",
                                                stringr::str_detect(geneSymbol,"SNORD116") ~ "SNORD116",
                                                TRUE ~ as.character(geneSymbol))) %>% 
    dplyr::select(geneSymbol) %>%
    dplyr::distinct()
  
  sigRegionsOverlap <- sigRegions %>%
    imprintOverlaps() %>%
    cleanOverlaps()
  
  regionsOverlap <- regions %>%
    imprintOverlaps() %>%
    cleanOverlaps()
  
  sigRegionsCount <- sigRegions %>%
    cleanOverlaps()
  
  regionsCount <- regions %>%
    cleanOverlaps()
  
  imprintMatrix <- matrix(c(nrow(sigRegionsOverlap), (nrow(sigRegionsCount) - nrow(sigRegionsOverlap)),
                            nrow(regionsOverlap), (nrow(regionsCount) - nrow(regionsOverlap))),
                          nrow = 2)
  
  p <- fisher.test(imprintMatrix, alternative = "greater")$p.value
  
  print(glue::glue("{nrow(sigRegionsOverlap)} out of {length(imprint)} human imprinted genes \\
                   are present in the {nrow(sigRegions)} DMRs
                   The p-value for the over-enrichment analysis is {round(p, digits = 2)}
                   The imprinted genes in the DMRs are: {tidyOverlaps}",
                   tidyOverlaps = glue::glue_collapse({dplyr::pull(sigRegionsOverlap)}, sep = ", ", last = " and ")))
                   
  
  return(cat("\n"))
}
