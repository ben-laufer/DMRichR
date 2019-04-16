#' GOplot
#' @description Plots top signficant Gene Ontology and KEGG pathway terms from enrichR
#' @param GO A list of ontology and pathway data frames returned from \code{enrichR::enrichr()}
#' @return A \code{ggplot} object of top significant GO and KEGG pathway terms from an enrichR analysis
#'  that can be viewed by calling it, saved with \code{ggplot2::ggsave()}, or further modified by adding \code{ggplot2} syntax
#' @import tidyverse
#' @import cowplot
#' @export GOplot

GOplot <- function(GO = GO){
  GOplot <- rbind(GO$GO_Biological_Process_2018[c(1:5),],
                  GO$GO_Cellular_Component_2018[c(1:5),],
                  GO$GO_Molecular_Function_2018[c(1:5),],
                  GO$KEGG_2016[c(1:5),]) %>%
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
    dplyr::mutate(Term = stringr::str_replace(.$Term, "_.*", "")) %>%
    dplyr::mutate(Term = stringr::str_trim(.$Term)) %>%
    dplyr::mutate(Term = stringr::str_to_title(.$Term)) %>%
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