#' imprintOverlap
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
#' @export imprintOverlap

imprintOverlap <- function(sigRegions = sigRegions,
                           regions = regions,
                           TxDb = TxDb,
                           annoDb = annoDb){
  
  imprint <- c("^DIRAS3$",
               "^RNU5D-1$",
               "^TP73$",
               "^LRRTM1$",
               "^GPR1$",
               "^ZDBF2$",
               "^NAP1L5$",
               "^ERAP2$",
               "^VTRNA2-1$",
               "^TNDM$",
               "^ADTRP$",
               "^FAM50B$",
               "^LIN28B$",
               "^AIM1$",
               "^PLAGL1$",
               "^HYMAI$",
               "^SLC22A2$",
               "^SLC22A3$",
               "^GRB10$",
               "^DDC$",
               "^PEG10$",
               "^MAGI2$",
               "^SGCE$",
               "^PPP1R9A$",
               "^TFPI2$",
               "^DLX5$",
               "^MEST$",
               "^CPA4$",
               "^COPG2IT1$",
               "^MESTIT1$",
               "^KLF14$",
               "^DLGAP2$",
               "^KCNK9$",
               "^ZFAT-AS1$",
               "^PEG13$",
               "^ZFAT$",
               "^GLIS3$",
               "^INPP5F$",
               "^WT1-AS$",
               "^WT1$",
               "^KCNQ1OT1$",
               "^OSBPL5$",
               "^KCNQ1DN$",
               "^KCNQ1$",
               "^H19$",
               "^CDKN1C$",
               "^PHLDA2$",
               "^IGF2AS$",
               "^SLC22A18$",
               "^INS$",
               "^IGF2$",
               "^ANO1$",
               "^ZC3H12C$",
               "^NTM$",
               "^RBP5$",
               "^LRP1$",
               "^ATP5F1EP2$",
               "^RB1$",
               "^SMOC1$",
               "^DLK1$",
               "^MEG3$",
               "^RTL1$",
               "^SNORD114",
               "^MEG8$",
               "^SNORD113",
               "^MAGEL2$",
               "^MKRN3$",
               "^UBE3A$",
               "^NPAP1$",
               "^SNORD109A$",
               "^SNORD108$",
               "^SNORD107$",
               "^SNORD109B$",
               "^ATP10A$",
               "^SNRPN$",
               "^SNORD116",
               "^SNORD115",
               "^PWCR1$",
               "^NDN$",
               "^SNURF$",
               "^SNORD64$",
               "^IRAIN$",
               "^ZNF597$",
               "^NAA60$",
               "^TP53$",
               "^TCEB3C$",
               "^DNMT1$",
               "^ZIM2$",
               "^PEG3$",
               "^MIMT1$",
               "^NLRP2$",
               "^MIR371A$",
               "^PSIMCT-1$",
               "^BLCAP$",
               "^NNAT$",
               "^MCTS2$",
               "^GDAP1L1$",
               "^SGK2$",
               "^GNAS$",
               "^L3MBTL$",
               "^SANG$",
               "^GNASAS$",
               "^MIR298$",
               "^MIR296$",
               "^DGCR6$",
               "^DGCR6L$")
  
  sigRegions <- sigRegions %>%
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
                                                TRUE ~ as.character(geneSymbol)
    )
    ) %>% 
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
                   The imprinted genes in the DMRs are:"))
                   print(glue::glue_collapse({dplyr::pull(sigRegionsOverlap)}, sep = ", "))
  
  return(cat("\n"))
}
