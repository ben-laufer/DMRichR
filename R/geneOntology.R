#' GOfuncR
#' @description Perform Gene Ontology enrichment analysis of DMRs using \code{GOfuncR}. 
#' @param sigRegions \code{GRanges} object of DMRs.
#' @param regions \code{GRanges} object of background regions. 
#' @param n_randsets Number specifying the number of random sets for calculating the FWER.
#' @param upstream Numeric of how many bases to extend upstream from gene body for mapping DMRs to genes.
#' @param downstream Numeric of how many bases to extend downstream from gene body for mapping DMRs to genes.
#' @param annoDb Character specifying \code{OrgDb} annotation package for species of interest.
#' @param TxDb \code{TxDb} or \code{EnsDb} annotation package for genome of interest.
#' @param ... Additional arguments passed onto \code{GOfuncR::go_enrich()}.
#' @import GOfuncR
#' @import GenomicRanges
#' @import ensembldb
#' @importFrom GenomicFeatures genes
#' @importFrom GenomeInfoDb keepStandardChromosomes as.data.frame seqlevelsStyle
#' @importFrom glue glue
#' @importFrom magrittr %>%
#' @importFrom dplyr as_tibble mutate distinct select
#' @importFrom tidyr unite
#' @importFrom plyranges count_overlaps
#' @export GOfuncR
#' 
GOfuncR <- function(sigRegions = sigRegions,
                    regions = regions,
                    n_randsets = 1000,
                    upstream = 5000,
                    downstream = 1000,
                    annoDb = annoDb,
                    TxDb = TxDb,
                    ...){
  
  print(glue::glue("Selecting annotation databases..."))
  
  if(is(TxDb, "TxDb")){
    
    print(glue::glue("Obtaining UCSC gene annotations..."))
    
    gene_coords <- TxDb %>%
      GenomicFeatures::genes() %>%
      GenomeInfoDb::keepStandardChromosomes(pruning.mode = "coarse") %>%
      DMRichR::extend(upstream = upstream, downstream = downstream) %>%
      dplyr::as_tibble() %>% 
      dplyr::mutate(gene_id = as.integer(.$gene_id)) %>% 
      dplyr::mutate(gene_id = GOfuncR:::entrez_to_symbol(.$gene_id, get(annoDb))[,2]) %>% 
      dplyr::distinct(gene_id, .keep_all = T) %>% 
      dplyr::select(symbol = gene_id, seqnames, start, end) %>% 
      as.data.frame() 
    
    print(glue::glue("{nrow(gene_coords)} unique genes will be utilized for GOfuncR..."))
    
  }else if(is(TxDb, "EnsDb")){
    
    print(glue::glue("Obtaining ENSEMBL gene annotations..."))
    
    genes <- TxDb %>%
      ensembldb::genes(., filter = GeneBiotypeFilter("protein_coding")) %>%
      GenomeInfoDb::keepStandardChromosomes(pruning.mode = "coarse")
    
    GenomeInfoDb::seqlevelsStyle(genes) <- "UCSC"
    
    gene_coords <- genes %>% 
      DMRichR::extend(upstream = upstream, downstream = downstream) %>% 
      dplyr::as_tibble() %>% 
      dplyr::distinct(gene_name, .keep_all = T) %>% 
      dplyr::select(symbol = gene_name, seqnames, start, end) %>% 
      dplyr::filter(symbol != "") %>% 
      as.data.frame()
    
    print(glue::glue("{nrow(gene_coords)} unique genes will be utilized for GOfuncR..."))
  }
  
  coord <- regions %>%
    plyranges::mutate(candidate = plyranges::count_overlaps(., sigRegions)) %>% 
    plyranges::mutate(candidate = dplyr::case_when(candidate != 0 ~ 1,
                                                   candidate == 0 ~ 0)) %>%
    GenomeInfoDb::as.data.frame() %>%
    dplyr::select(seqnames, start, end, candidate) %>%
    tidyr::unite(c("seqnames","start"), col = "seqstart", sep = ":") %>%
    tidyr::unite(c("seqstart","end"), col = "coordinate", sep = "-") 
  
  coord$candidate %>%
    table()
  
  print(glue::glue("Performing enrichment testing..."))
  
  GOfuncResults <- GOfuncR::go_enrich(genes = coord,
                                      test = 'hyper',
                                      n_randsets = n_randsets,
                                      regions = TRUE,
                                      gene_coords = gene_coords,
                                      circ_chrom = TRUE, # Otherwise get the error: "Background regions too small."
                                      gene_len = TRUE,
                                      orgDb = annoDb, # Blocking this makes it use human mappings, which are better in most cases for non-model organism since their org.dbs contain old GO mappings
                                      #txDb = TxDb, # Not used for custom gene coords
                                      silent = TRUE,
                                      ...)
  
  return(GOfuncResults)
}

#' slimGO
#' @description Slims top significant Gene Ontology terms from \code{enrichR}, \code{rGREAT}, 
#' and \code{GOfuncR} using \code{rrvgo}.
#' @param GO A dataframe or list of dataframes returned from \code{enrichR::enrichr()}, 
#' \code{rGREAT::getEnrichmentTables()},or \code{GOfuncR::go_enrich()}.
#' @param tool A character vector of the name of the database (enrichR, rGREAT, or GOfuncR).
#' @param annoDb Character specifying \code{OrgDb} annotation package for species of interest.
#' @param plots Logical indicating if scatter and treemap plots should be generated.
#' @return A \code{tibble} of top distinct and significant GO terms from an \code{enrichR},
#'  \code{rGREAT} or \code{GOfuncR} analysis.
#' @import enrichR
#' @import rGREAT
#' @import GOfuncR
#' @import rrvgo
#' @importFrom magrittr %>%
#' @importFrom dplyr filter as_tibble mutate select
#' @importFrom data.table rbindlist
#' @importFrom glue glue
#' @importFrom purrr set_names map_dfr
#' @export slimGO
#' 
slimGO <- function(GO = GO,
                   tool = c("enrichR", "rGREAT", "GOfuncR"),
                   annoDb = annoDb,
                   plots = FALSE){
  
  if(tool == "enrichR"){
    GO <- GO %>%
      data.table::rbindlist(idcol = "Gene Ontology") %>%
      dplyr::as_tibble() %>%
      dplyr::filter(`Gene Ontology` %in% c("GO_Biological_Process_2018",
                                           "GO_Cellular_Component_2018",
                                           "GO_Molecular_Function_2018")) %>% 
      dplyr::mutate(Term = stringr::str_extract(.$Term, "\\(GO.*")) %>%
      dplyr::mutate(Term = stringr::str_replace_all(.$Term, "[//(//)]",""), "") %>%
      dplyr::mutate("Gene Ontology" = dplyr::case_when(`Gene Ontology` == "GO_Biological_Process_2018" ~ "BP",
                                                       `Gene Ontology` == "GO_Cellular_Component_2018" ~ "CC",
                                                       `Gene Ontology` == "GO_Molecular_Function_2018" ~ "MF")) %>%
      dplyr::select(p = P.value, go = Term, "Gene Ontology") %>% 
      dplyr::filter(p < 0.05)
    
  }else if(tool == "rGREAT"){
    
    GO <-  GO %>%
      data.table::rbindlist(idcol = "Gene Ontology") %>%
      dplyr::as_tibble() %>%
      dplyr::mutate("Gene Ontology" = dplyr::case_when(`Gene Ontology` == "GO Biological Process" ~ "BP",
                                                       `Gene Ontology` == "GO Cellular Component" ~ "CC",
                                                       `Gene Ontology` == "GO Molecular Function" ~ "MF")) %>%
      dplyr::select(p = Hyper_Raw_PValue, go = ID, "Gene Ontology") %>% 
      dplyr::filter(p < 0.05)
    
  }else if(tool == "GOfuncR"){
    
    GO <- GO$results %>%
      dplyr::as_tibble() %>%
      dplyr::mutate("Gene Ontology" = dplyr::case_when(ontology == "biological_process" ~ "BP",
                                                       ontology == "cellular_component" ~ "CC",
                                                       ontology == "molecular_function" ~ "MF")) %>%
      dplyr::select(p = raw_p_overrep, go = node_id, "Gene Ontology") %>% 
      dplyr::filter(p < 0.05)
    
  }else{
    stop(glue("{tool} is not supported, please choose either enrichR, rGREAT, or GOfuncR [Case Sensitive]"))
  }
  
  print(glue::glue("Submiting results from {tool} to rrvgo..."))
  
  .slim <- function(GO = GO,
                    ont = ont,
                    annoDb = annoDb,
                    plots = plots,
                    tool = tool,
                    threshold = threshold){
    GO <- GO %>%
      dplyr::filter(`Gene Ontology` == ont)
    
    print(glue::glue("rrvgo is now slimming {ont} GO terms from {tool}"))
    
    simMatrix  <- rrvgo::calculateSimMatrix(GO$go,
                                            orgdb = annoDb,
                                            ont = ont,
                                            method = "Rel")
    
    reducedTerms <- rrvgo::reduceSimMatrix(simMatrix,
                                           setNames(-log10(GO$p), GO$go),
                                           threshold = threshold,
                                           orgdb = annoDb) 
    
    if(plots == TRUE){
      p <- rrvgo::scatterPlot(simMatrix, reducedTerms) # Doesn't plot otherwise
      plot(p) 
      rrvgo::treemapPlot(reducedTerms)
    }
    
    print(glue::glue("There are {max(reducedTerms$cluster)} clusters in your GO {ont} terms from {tool}"))
    
    reducedTerms %>%   
      dplyr::as_tibble() %>%
      return()
  }
  
  slimmed <- GO %>%
    dplyr::select(`Gene Ontology`) %>%
    table() %>%
    names() %>% 
    purrr::set_names() %>%
    purrr::map_dfr(~.slim(GO = GO,
                          ont = .,
                          annoDb = annoDb,
                          tool = tool,
                          plots = plots,
                          threshold = threshold),
                   .id = "Gene Ontology") %>%
    dplyr::inner_join(GO) %>%
    dplyr::filter(term == as.character(parentTerm)) %>%
    dplyr::mutate("-log10.p-value" = -log10(p)) %>%
    dplyr::mutate("Gene Ontology" = dplyr::recode_factor(`Gene Ontology`,
                                                         "BP" = "Biological Process",
                                                         "CC" = "Cellular Component",
                                                         "MF" = "Molecular Function")) %>%
    dplyr::arrange(dplyr::desc(`-log10.p-value`)) %>% 
    dplyr::select("Gene Ontology", Term = term, "-log10.p-value") %>% 
    return()
}

#' GOplot
#' @description Plots top significant slimmed Gene Ontology terms from from \code{enrichR},
#'  \code{rGREAT}, and \code{GOfuncR}.
#' @param slimmedGO A \code{tibble} from \code{DMRichR::rrvgo()} or \code{DMRichR::REVIGO()}.
#' @return A \code{ggplot} object of top significant GO and pathway terms from an \code{enrichR}, 
#' \code{GOfuncR}, or \code{rGREAT} analysis that can be viewed by calling it, saved with 
#' \code{ggplot2::ggsave()}, or further modified by adding \code{ggplot2} syntax.
#' @import ggplot2
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate select group_by slice ungroup
#' @importFrom forcats fct_rev
#' @importFrom ggsci scale_fill_d3
#' @importFrom Hmisc capitalize
#' @importFrom stringr str_trim str_trunc
#' @importFrom glue glue
#' @export GOplot
#' 
GOplot <- function(slimmedGO = slimmedGO){
  
  slimmedGO %>% 
    dplyr::group_by(`Gene Ontology`) %>%
    dplyr::slice(1:7) %>%
    dplyr::ungroup() %>% 
    dplyr::mutate(Term = stringr::str_trim(.$Term)) %>%
    #dplyr::mutate(Term = Hmisc::capitalize(.$Term)) %>%
    dplyr::mutate(Term = stringr::str_trunc(.$Term, 40, side = "right")) %>% 
    dplyr::mutate(Term = factor(.$Term, levels = unique(.$Term[order(forcats::fct_rev(.$`Gene Ontology`), .$`-log10.p-value`)]))) %>% 
    ggplot2::ggplot(aes(x = Term,
                        y = `-log10.p-value`,
                        fill = `Gene Ontology`,
                        group = `Gene Ontology`)) +
    ggplot2::geom_bar(stat = "identity",
                      position = position_dodge(),
                      color = "Black") +
    ggplot2::coord_flip() +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggsci::scale_fill_d3() +
    ggplot2::labs(y = expression("-log"[10](p))) +
    ggplot2::theme_classic() +
    ggplot2::theme(text = element_text(size = 40),
                   axis.title.y = element_blank(),
                   axis.title.x = element_text(size = 25)) %>% 
    return()
}
