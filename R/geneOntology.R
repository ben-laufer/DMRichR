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
#' @description Slims and plots top signficant Gene Ontology terms from enrichR, rGREAT, and GOfuncR.
#' The terms are ranked by dispensability before being plotted, which then orders them by p-value. 
#' @param GO A dataframe or list of dataframes returned
#' from \code{enrichR::enrichr()}, \code{rGREAT::getEnrichmentTables()}, or \code{GOfuncR::go_enrich()}.
#' @param tool A character vector of the name of the database (enrichR, rGREAT, or GOfuncR).
#' @return A \code{ggplot} object of top significant GO and pathway terms from an \code{enrichR} 
#' or \code{rGREAT} analysis that can be viewed by calling it, saved with \code{ggplot2::ggsave()}, 
#' or further modified by adding \code{ggplot2} syntax.
#' @import tidyverse
#' @import enrichR
#' @import rGREAT
#' @import GOfuncR
#' @import ggplot2
#' @importFrom ggsci scale_fill_d3
#' @importFrom rvest html_session html_form set_values submit_form
#' @importFrom Hmisc capitalize
#' @importFrom data.table rbindlist
#' @importFrom glue glue
#' @references \url{https://github.com/hbc/revigoR/blob/master/rvest_revigo.R}
#' @references \url{http://revigo.irb.hr}
#' @export GOplot

GOplot <- function(GO = GO,
                   tool = c("enrichR", "rGREAT", "GOfuncR")){
  
  print(glue::glue("Tidying results from {tool}..."))
  if(tool == "enrichR"){
    GO <- GO %>%
      data.table::rbindlist(idcol = "Database") %>%
      dplyr::filter(Database %in% c("GO_Biological_Process_2018",
                                    "GO_Cellular_Component_2018",
                                    "GO_Molecular_Function_2018"
                                    )
                    ) %>% 
      dplyr::as_tibble() %>%
      dplyr::mutate(Term = stringr::str_extract(.$Term, "\\(GO.*")) %>%
      dplyr::mutate(Term = stringr::str_replace_all(.$Term, "[//(//)]",""), "") %>%
      dplyr::filter(P.value <= 0.05)
    
    goList <- paste(GO$Term, GO$P.value, collapse = "\n")
    
  }else if(tool == "rGREAT"){
    
    GO <-  GO %>%
      data.table::rbindlist(idcol = "Database") %>%
      dplyr::as_tibble() %>%
      dplyr::filter(Hyper_Raw_PValue <= 0.05)
    
    goList <- paste(GO$ID, GO$Hyper_Raw_PValue, collapse = "\n")
    
  }else if(tool == "GOfuncR"){
    
    GO <- GO$results %>%
      dplyr::filter(raw_p_overrep <= 0.05)
    
    goList <- paste(GO$node_id, GO$raw_p_overrep, collapse = "\n")
    
  }else{
    stop(glue("{tool} is not supported, please choose either enrichR, rGREAT, or GOfuncR [Case Sensitive]"))
  }
  
  print(glue::glue("Submiting results from {tool} to REVIGO..."))
  revigo_session <- rvest::html_session("http://revigo.irb.hr/")
  revigo_form <- rvest::html_form(revigo_session)[[1]]  
  filled_form <- rvest::set_values(revigo_form,
                                   'goList' = goList,
                                   'cutoff' = 0.4,
                                   'isPValue' = "yes",
                                   'measure' = "SIMREL")
  result_page <- rvest::submit_form(revigo_session,
                                    filled_form,
                                    submit = 'startRevigo')
  
  revigo_results <- list()
  for (i in 1:3){
    results_table <- rvest::html_table(result_page)[[i]]
    names(results_table) <- results_table[2,]
    revigo_results[[i]] <- results_table[3:nrow(results_table),]
  }
  names(revigo_results) <- c("Biological Process", "Cellular Component", "Molecular Function")
  
  revigo_results <- data.table::rbindlist(revigo_results, idcol = "Gene Ontology") %>%
    dplyr::filter(dispensability < 0.4) %>%
    dplyr::select("Gene Ontology",
                  Term = "description",
                  "-log10.p-value" = `log10 p-value`) %>% 
    dplyr::mutate("-log10.p-value" = -(as.numeric(`-log10.p-value`))) %>% 
    #dplyr::arrange(dplyr::desc(`-log10.p-value`)) %>% 
    dplyr::mutate("Gene Ontology" = as.factor(`Gene Ontology`)) %>% 
    dplyr::group_by(`Gene Ontology`) %>%
    dplyr::slice(1:7) %>%
    dplyr::ungroup() %>% 
    dplyr::mutate(Term = stringr::str_trim(.$Term)) %>%
    dplyr::mutate(Term = Hmisc::capitalize(.$Term)) %>%
    dplyr::mutate(Term = stringr::str_wrap(.$Term, 45)) %>% 
    dplyr::mutate(Term = factor(.$Term, levels = unique(.$Term[order(forcats::fct_rev(.$`Gene Ontology`), .$`-log10.p-value`)]))) 

  print(glue::glue("Plotting slimmed gene ontology results from {tool}..."))
  GOplot <- revigo_results %>%
    ggplot2::ggplot(aes(x = Term, y = `-log10.p-value`, fill = `Gene Ontology`, group = `Gene Ontology`)) +
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