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
#' @export GOplot

GOplot <- function(GO = GO,
                   tool = NULL){
  
  if(tool == "enrichR"){
    GOplot <- rbind(GO$GO_Biological_Process_2018[c(1:5),],
                    GO$GO_Cellular_Component_2018[c(1:5),],
                    GO$GO_Molecular_Function_2018[c(1:5),],
                    GO$KEGG_2016[c(1:5),]
                    ) %>%
      dplyr::as_tibble() %>%
      cbind(
        dplyr::as_tibble(
          c(
            rep("Biological Process", 5),
            rep("Cellular Component", 5),
            rep("Molecular Function", 5),
            rep("KEGG", 5)
          )
        )
      ) %>%
      dplyr::select(Term, P.value, value, Combined.Score) %>%
      dplyr::filter(P.value <= 0.05) %>%
      dplyr::mutate(P.value = -log10(P.value)) %>%
      dplyr::rename(`-log10.p-value` = P.value) %>%
      dplyr::rename(Database = value) %>% 
      dplyr::mutate(Database = stringr::str_replace(.$Database, "KEGG", "Pathway (KEGG)")) %>% 
      dplyr::mutate(Term = stringr::str_replace(.$Term, "\\(.*", "")) %>%
      dplyr::mutate(Term = stringr::str_replace(.$Term, "_.*", ""))
  
  }else if(tool == "rGREAT"){
    GOplot <- rbind(GO$`GO Biological Process`[c(1:5),],
                    GO$`GO Cellular Component`[c(1:5),],
                    GO$`GO Molecular Function`[c(1:5),],
                    GO$`MSigDB Pathway`[c(1:5),]
                    ) %>%
      dplyr::as_tibble() %>%
      cbind(
        dplyr::as_tibble(
          c(
            rep("Biological Process", 5),
            rep("Cellular Component", 5),
            rep("Molecular Function", 5),
            rep("MSigDB Pathway", 5)
          )
        )
      ) %>%
      dplyr::select(name, Hyper_Raw_PValue, value) %>%
      dplyr::filter(Hyper_Raw_PValue <= 0.05) %>%
      dplyr::mutate(Hyper_Raw_PValue = -log10(Hyper_Raw_PValue)) %>%
      dplyr::rename(`-log10.p-value`= Hyper_Raw_PValue,
                    Database = value,
                    Term = name) 
  
  }else if(tool == "GOfuncR"){
    GOplot <- rbind(GO$results %>% dplyr::filter(ontology == "biological_process") %>% dplyr::slice(1:5),
                    GO$results %>% dplyr::filter(ontology == "cellular_component") %>% dplyr::slice(1:5),
                    GO$results %>% dplyr::filter(ontology == "molecular_function") %>% dplyr::slice(1:5)
                    ) %>%
      dplyr::as_tibble() %>%
      cbind(
        dplyr::as_tibble(
          c(
            rep("Biological Process", 5),
            rep("Cellular Component", 5),
            rep("Molecular Function", 5)
          )
        )
      ) %>%
      dplyr::select(node_name, raw_p_overrep, value) %>%
      dplyr::filter(raw_p_overrep <= 0.05) %>%
      dplyr::mutate(raw_p_overrep = -log10(raw_p_overrep)) %>%
      dplyr::rename(`-log10.p-value`= raw_p_overrep,
                    Database = value,
                    Term = node_name) 
   
  }else{
    stop(glue("{tool} is not supported, please choose either enrichR or rGREAT [Case Sensitive]"))
  }
   
  GOplot <- GOplot %>%
    dplyr::mutate(Term = stringr::str_trim(.$Term)) %>%
    dplyr::mutate(Term = stringr::str_to_title(.$Term)) %>%
    dplyr::mutate(Term = stringr::str_wrap(.$Term, 70)) %>% 
    dplyr::mutate(Database = factor(.$Database)) %>% 
    dplyr::mutate(Term = factor(.$Term, levels = unique(.$Term[order(forcats::fct_rev(.$Database), .$`-log10.p-value`)]))) %>% 
    ggplot2::ggplot(aes(x = Term, y = `-log10.p-value`, fill = Database, group = Database)) +
    geom_bar(stat = "identity", position = position_dodge(), color = "Black") +
    coord_flip() +
    scale_y_continuous(expand = c(0, 0)) +
    labs(y = expression("-log"[10](p))) +
    theme_classic() +
    theme(axis.text = element_text(size = 14),
          axis.title.y = element_blank())
  
  return(GOplot)
}