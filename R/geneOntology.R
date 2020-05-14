#' GOfuncR
#' @description Perform Gene Ontology enrichment analysis of DMRs using \code{GOfuncR}. 
#' @param sigRegions \code{GRanges} object of DMRs.
#' @param regions \code{GRanges} object of background regions. 
#' @param n_randsets Number specifying the number of random sets for calculating the FWER.
#' @param upstream Numeric of how many bases to extend upstream from gene body for mapping DMRs to genes.
#' @param downstream Numeric of how many bases to extend downstream from gene body for mapping DMRs to genes.
#' @param annoDb Character specifying OrgDb annotation package for species of interest.
#' @param TxDb TxDb annotation package for genome of interest.
#' @param ... Additional arugments passed onto \code{GOfuncR::go_enrich()}.
#' @import GOfuncR
#' @import tidyverse
#' @importFrom glue glue
#' @references \url{https://support.bioconductor.org/p/78652/}
#' @export GOfuncR
GOfuncR <- function(sigRegions = sigRegions,
                    regions = regions,
                    n_randsets = 1000,
                    upstream = 5000,
                    downstream = 1000,
                    annoDb = annoDb,
                    TxDb = TxDb,
                    ...){
  
  cat("\n[DMRichR] GOfuncR \t\t\t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
  
  print(glue::glue("Selecting annotation databases..."))

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
    dplyr::mutate(gene_id = GOfuncR:::entrez_to_symbol(.$gene_id, get(annoDb))[,2]) %>% 
    dplyr::distinct(gene_id, .keep_all = T) %>% 
    dplyr::select(symbol = gene_id, seqnames, start, end) %>% 
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
                                      silent = TRUE,
                                      ...)
  
  return(GOfuncResults)
  
}

#' GOplot
#' @description Plots top signficant Gene Ontology and pathway terms from enrichR, rGREAT, and GOfuncR.
#' @param GO A list of ontology and pathway data frames returned from \code{enrichR::enrichr()}, \code{rGREAT::getEnrichmentTables()}, or \code{GOfuncR::go_enrich()}.
#' @param tool A character vector of the name of the database (enrichR, rGREAT, or GOfuncR).
#' @return A \code{ggplot} object of top significant GO and pathway terms from an \code{enrichR} or \code{rGREAT} analysis.
#'  that can be viewed by calling it, saved with \code{ggplot2::ggsave()}, or further modified by adding \code{ggplot2} syntax.
#' @import tidyverse
#' @import enrichR
#' @import rGREAT
#' @import GOfuncR
#' @import ggplot2
#' @import ggsci
#' @importFrom Hmisc capitalize
#' @importFrom data.table rbindlist
#' @export GOplot

GOplot <- function(GO = GO,
                   tool = NULL){
  
  if(tool == "enrichR"){
   GOplot <- GO %>%
     data.table::rbindlist(idcol = "Database") %>%
     dplyr::filter(Database %in% c("GO_Biological_Process_2018",
                                   "GO_Cellular_Component_2018",
                                   "GO_Molecular_Function_2018",
                                   "KEGG_2019_Human"
                                   )
                   ) %>% 
     dplyr::as_tibble() %>%
     dplyr::mutate(Database = dplyr::recode_factor(Database,
                                                   "GO_Biological_Process_2018" = "Biological Process",
                                                   "GO_Cellular_Component_2018" = "Cellular Component",
                                                   "GO_Molecular_Function_2018" = "Molecular Function",
                                                   "KEGG_2019_Human" = "KEGG (Pathway)"
                                                   )
                   ) %>% 
     dplyr::select(Term, P.value, Database) %>%
     dplyr::filter(P.value <= 0.05) %>%
     dplyr::mutate(P.value = -log10(P.value)) %>%
     dplyr::rename(`-log10.p-value` = P.value) %>%
     dplyr::mutate(Term = stringr::str_replace(.$Term, "\\(GO.*", "")) %>%
     dplyr::mutate(Term = stringr::str_replace(.$Term, "_.*", "")) %>% 
     dplyr::group_by(Database) %>%
     dplyr::slice(1:7) %>%
     dplyr::ungroup()
    
  }else if(tool == "rGREAT"){
    GOplot <- GO %>%
      data.table::rbindlist(idcol = "Database") %>%
      dplyr::as_tibble() %>%
      dplyr::mutate(Database = dplyr::recode_factor(Database,
                                                    "GO Biological Process" = "Biological Process",
                                                    "GO Cellular Component" = "Cellular Component",
                                                    "GO Molecular Function" = "Molecular Function"
                                                    )
                    ) %>% 
      dplyr::filter(Hyper_Raw_PValue <= 0.05) %>%
      dplyr::mutate(Hyper_Raw_PValue = -log10(Hyper_Raw_PValue)) %>%
      dplyr::select(Term = "name", `-log10.p-value` = Hyper_Raw_PValue, Database) %>%
      dplyr::group_by(Database) %>%
      dplyr::slice(1:7) %>%
      dplyr::ungroup()
    
  }else if(tool == "GOfuncR"){
    GOplot <- GO$results %>%
      dplyr::mutate(ontology = dplyr::recode_factor(ontology,
                                                    "biological_process" = "Biological Process",
                                                    "cellular_component" = "Cellular Component",
                                                    "molecular_function" = "Molecular Function")
                    ) %>%
      dplyr::filter(raw_p_overrep <= 0.05) %>%
      dplyr::mutate(raw_p_overrep = -log10(raw_p_overrep)) %>%
      dplyr::group_by(ontology) %>%
      dplyr::slice(1:7) %>%
      dplyr::ungroup() %>% 
      dplyr::select(Term = node_name,
                    Database = ontology,
                    "-log10.p-value"= raw_p_overrep)
   
  }else{
    stop(glue("{tool} is not supported, please choose either enrichR, rGREAT, or GOfuncR [Case Sensitive]"))
  }
   
  GOplot <- GOplot %>%
    dplyr::mutate(Term = stringr::str_trim(.$Term)) %>%
    dplyr::mutate(Term = Hmisc::capitalize(.$Term)) %>%
    dplyr::mutate(Term = stringr::str_wrap(.$Term, 45)) %>% 
    dplyr::mutate(Database = factor(.$Database)) %>% 
    dplyr::mutate(Term = factor(.$Term, levels = unique(.$Term[order(forcats::fct_rev(.$Database), .$`-log10.p-value`)]))) %>% 
    ggplot2::ggplot(aes(x = Term, y = `-log10.p-value`, fill = Database, group = Database)) +
    geom_bar(stat = "identity", position = position_dodge(), color = "Black") +
    coord_flip() +
    scale_y_continuous(expand = c(0, 0)) +
    ggsci::scale_fill_d3() +
    labs(y = expression("-log"[10](p))) +
    theme_classic() +
    theme(axis.text = element_text(size = 14),
          axis.title.y = element_blank(),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 14),
          strip.text = element_text(size = 14)
    ) %>% 
    return()
}