#' GOfuncR
#' @description Perform Gene Ontology enrichment analysis of DMRs using \code{GOfuncR}. 
#' @param sigRegions \code{GRanges} object of DMRs.
#' @param regions \code{GRanges} object of background regions. 
#' @param genome Character specifying genome of interest ("hg38", "mm10", rheMac8", "rn6").
#' @param n_randsets Number specifying the number of random sets for calculating the FWER.
#' @param upstream Numeric of how many bases to extend upstream from gene body for mapping DMRs to genes.
#' @param downstream Numeric of how many bases to extend downstream from gene body for mapping DMRs to genes.
#' @param annoDb Character specifying OrgDb annotation package for species of interest.
#' @param TxDb TxDb annotation package for genome of interest.
#' @param ... Additional arugments passed onto \code{GOfuncR::go_enrich()}.
#' @import tidyverse
#' @import annotables
#' @export GOfuncR
GOfuncR <- function(sigRegions = sigRegions,
                    regions = regions,
                    genome = NULL,
                    n_randsets = 1000,
                    upstream = 5000,
                    downstream = 1000,
                    annoDb = annoDb,
                    TxDb = TxDb,
                    ...){
  
  cat("\n[DMRichR] GOfuncR \t\t\t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
  
  print(glue::glue("Selecting annotation databases..."))
  # Error: Can't use numeric NA as column index with `[`.
  # annotable <- dplyr::case_when(genome == "hg38" ~ annotables::grch38,
  #                               genome == "mm10" ~ annotables::grcm38,
  #                               genome == "rheMac8" ~ annotables::mmul801,
  #                               genome == "rn6" ~ annotables::rnor6
  # )
  
  if(genome == "hg38"){
    annotable <- annotables::grch38
  }else if(genome == "mm10"){
    annotable <- annotables::grcm38
  }else if(genome == "rheMac8"){
    annotable <- annotables::mmul801
  }else if(genome == "rn6"){
    annotable <- annotables::rnor6
  }else{
    stop(glue("{genome} is not supported, please choose either hg38, mm10, rheMac8, or rn6 [Case Sensitive]"))
  }
  
  #https://support.bioconductor.org/p/78652/
  extend <- function(x,
                     upstream = 0,
                     downstream = 0)
  {
    if (any(strand(x) == "*"))
      warning("'*' ranges were treated as '+'")
    on_plus <- strand(x) == "+" | strand(x) == "*"
    new_start <- start(x) - ifelse(on_plus, upstream, downstream)
    new_end <- end(x) + ifelse(on_plus, downstream, upstream)
    ranges(x) <- IRanges(new_start, new_end)
    trim(x)
  }
  
  gene_coords <- TxDb %>%
    GenomicFeatures::genes() %>%
    GenomeInfoDb::keepStandardChromosomes(pruning.mode = "coarse") %>%
    extend(upstream = upstream, downstream = downstream) %>%
    dplyr::as_tibble() %>% 
    dplyr::mutate(gene_id = as.integer(.$gene_id)) %>% 
    dplyr::inner_join(annotable, by = c("gene_id" = "entrez")) %>%
    dplyr::distinct(.$symbol, .keep_all = T) %>%
    dplyr::select(symbol, seqnames, start.x, end.x) %>%
    dplyr::rename(start = start.x, end = end.x) %>%
    as.data.frame()
  
  ranges2coord <- . %>%
    GenomeInfoDb::as.data.frame() %>%
    dplyr::select(seqnames, start, end) %>%
    tidyr::unite(c("seqnames","start"), col = "seqstart", sep = ":") %>%
    tidyr::unite(c("seqstart","end"), col = "coordinate", sep = "-") %>%
    dplyr::as_tibble()
  
  sigRegions_coord <- sigRegions %>%
    ranges2coord()
  
  regions_coord <- regions %>%
    ranges2coord()
  
  coord <- data.frame(regions_coord,
                      candidate = 0)
  
  coord$candidate[which(sigRegions_coord$coordinate %in% regions_coord$coordinate)] <- 1
  
  coord$candidate %>%
    table()
  
  print(glue::glue("Performing enrichment testing..."))
  
  GOfuncResults <- GOfuncR::go_enrich(genes = coord,
                                      test = 'hyper',
                                      n_randsets = n_randsets,
                                      regions = TRUE,
                                      gene_coords = gene_coords,
                                      circ_chrom = TRUE, # Otherwise get the error: "Background regions too small."
                                      orgDb = annoDb,
                                      txDb = TxDb,
                                      silent = FALSE,
                                      ...)
  
  return(GOfuncResults)
  
}

# Alternate approach for annotations
# gene_coords <- annotable %>%
#   dplyr::filter(biotype == "protein_coding") %>%
#   dplyr::filter(chr %in% (regions %>%
#                             seqlevels() %>%
#                             stringr::str_remove("chr")
#   )
#   )%>% 
#   dplyr::distinct(.$symbol, .keep_all = T) %>%
#   dplyr::rename(gene = symbol, chromosome = chr) %>%
#   dplyr::select(gene, chromosome, start, end, strand) %>%
#   as.data.frame() %>%
#   GenomicRanges::makeGRangesFromDataFrame(.,
#                                           ignore.strand = FALSE,
#   )
