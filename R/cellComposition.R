#' bsseqLift
#' @description LiftOver a hg38 bsseq objet to hg19 coordinates
#' @param bsseq A \code{bsseq}object with hg38 coordinates
#' @return A \code{bsseq}object with hg19 coordinates
#' @import bsseq
#' @import tidyverse
#' @import rtracklayer
#' @import AnnotationHub
#' @import GenomicRanges
#' @importFrom glue glue
#' @export bsseqLift
bsseqLift <- function(bsseq = bs.filtered.bsseq){
  # Make indices
  mcols(bs.filtered.bsseq)$index <- 1:length(bs.filtered.bsseq)
  hg38 <- rowRanges(bs.filtered.bsseq)
  hg38$index <- 1:length(hg38)
  
  # Liftover ranges
  hg19 <- hg38 %>%
    rtracklayer::liftOver(AnnotationHub::AnnotationHub()[["AH14108"]]) %>%
    unlist()
  
  # Subset based on regions that lifted over
  bs.filtered.bsseq.hg19 <- bs.filtered.bsseq[which(hg19$index %in% hg38$index)]
  
  # Assign lifted over regions
  rowRanges(bs.filtered.bsseq.hg19) <- hg19
  
  print(glue::glue("{length(bs.filtered.bsseq.hg19)} out of {length(bs.filtered.bsseq)} were lifted over"))
  print(glue::glue("{length(bs.filtered.bsseq) - length(bs.filtered.bsseq.hg19)} did not liftOver"))
  
  genome(bs.filtered.bsseq.hg19) <- "hg19"
  return(bs.filtered.bsseq.hg19)
}

#' arrayRanges
#' @description Obtain hg19 EPIC array coordinates
#' @return A \code{GRanges} object of hg19 EPIC coordinates
#' @import tidyverse
#' @import GenomicRanges
#' @importFrom minfi getAnnotation
#' @references \url{https://support.bioconductor.org/p/78652/}
#' @export arrayRanges
arrayRanges <- function(){
  
  if(!require(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)){
    BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")}
  library("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
  
  message("Fetching coordinates for hg19...")
  
  array <- minfi::getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19) %>%
    as.data.frame() %>% 
    tibble::rownames_to_column() %>% 
    dplyr::select(rowname, chr, pos, strand) %>% 
    GenomicRanges::makeGRangesFromDataFrame(.,seqnames.field = "chr",
                                            start.field = "pos",
                                            end.field = "pos",
                                            strand.field = "strand",
                                            keep.extra.columns = TRUE) 
  
  names(array) <- array$rowname  
  array <- array %>%
    granges()
  
  return(array)
}

#' CCstats
#' @description Computes the cell composition differences between groups while adjusting for the provided covariates. 
#'  The differences are tested for using an ANOVA through the \code{\link[stats]{aov}} function.
#' @param CC A \code{matrix} from \code{FlowSorted.Blood.EPIC::projectCellType_CP()}
#' @param bsseq Smoothed bsseq object with design matrix in pData
#' @param testCovariate The factor to test for differences between groups
#' @param adjustCovariate The covariate(s) to adjust for between groups
#' @param matchCovariate Another covariate to adjust for between groups (for dmrseq compatibility)
#' @return A list of tibbles with the statsitics and the values used for the tests
#' @references \url{https://cran.r-project.org/web/packages/broom/vignettes/broom_and_dplyr.html}
#' @import bsseq
#' @import GenomeInfoDb
#' @import GenomicRanges
#' @import tidyverse
#' @import broom
#' @export CCstats
CCstats <- function(CC = CC,
                    bsseq = bs.filtered.bsseq,
                    testCovariate = testCovariate,
                    adjustCovariate = NULL,
                    matchCovariate = NULL){
  
  # Tidy --------------------------------------------------------------------
  
  tidyCC <- CC %>% 
    "*"(100) %>%
    "/"(rowSums(.)) %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column("Sample") %>% 
    dplyr::as_tibble() %>% 
    dplyr::full_join(bs.filtered.bsseq %>%
                       pData() %>%
                       as.data.frame() %>% 
                       tibble::rownames_to_column("Sample") %>%
                       dplyr::as_tibble(),
                     .)
  
  summary <- tidyCC %>%
    dplyr::select(one_of(!!testCovariate, !!adjustCovariate, !!matchCovariate, !!IDs)) %>% 
    dplyr::group_by_(testCovariate) %>%
    dplyr::summarise_at(IDs, mean)
  
  # Stats -------------------------------------------------------------------
  
  cat("Selecting model...")
  
  if(is.null(adjustCovariate) &
     (is.null(matchCovariate) | (length(levels(matchCovariate))) <= 1)){
    model <- as.formula(paste0("cellCount ~ ", paste(testCovariate)))
    
  }else if(!is.null(adjustCovariate) &
           (is.null(matchCovariate) | (length(levels(matchCovariate))) <= 1)){
    model <- as.formula(paste0("cellCount ~ ", paste(testCovariate, "+"), paste(adjustCovariate, collapse = " + ")))
    
  }else if(is.null(adjustCovariate) &
           (!is.null(matchCovariate) | !(length(levels(matchCovariate))) <= 1)){
    model <- as.formula(paste0("cellCount ~ ", paste(testCovariate, "+"), paste(matchCovariate)))
    
  }else if(!is.null(adjustCovariate) &
           (!is.null(matchCovariate) | !(length(levels(matchCovariate))) <= 1)){
    model <- as.formula(paste0("cellCount ~ ", paste(testCovariate, "+"), paste(adjustCovariate, collapse = " + "), paste(" + ", matchCovariate)))
  }
  cat("Done", "\n")
  cat(paste("The model is", paste(capture.output(print(model))[1], collapse= ' ')), "\n")
  
  IDs <- c("Neu", "NK", "Bcell" , "CD4T", "CD8T", "Mono")
  
  ANOVA <- tidyCC %>%
    tidyr::pivot_longer(cols = all_of(IDs),
                        names_to = "cellType",
                        values_to = "cellCount") %>%
    tidyr::nest(-cellType) %>%
    dplyr::mutate(
      ANOVA = purrr::map(data, ~ aov(model, data = .x)),
      tidied = purrr::map(ANOVA, broom::tidy)
    ) %>%
    dplyr::select(cellType, tidied) %>% 
    tidyr::unnest(tidied)
  
  list("input" = tidyCC,
       "summary" = summary,
       "ANOVA" = ANOVA) %>%
    return()
  
  pairWise <- tidyCC %>%
    tidyr::pivot_longer(cols = all_of(IDs),
                        names_to = "cellType",
                        values_to = "cellCount")
    tidyr::nest(-cellType) %>%
    dplyr::mutate(
      pairWise = purrr::map(data, ~ lm(model, data = .x) %>% 
                              lsmeans::ref.grid(data = .x) %>%
                              lsmeans::lsmeans(as.formula(paste("~", testCovariate))) %>% 
                              pairs(reverse = TRUE) %>%
                              summary()
      )
    ) %>%
    dplyr::select(cellType, pairWise) %>% 
    tidyr::unnest(pairWise) %>%
    dplyr::mutate(fdr = p.adjust(p.value, method = 'fdr'))
  
  list("input" = tidyCC,
       "summary" = summary,
       "ANOVA" = ANOVA,
       "pairWise" = pairWise
       ) %>%
    return()
  
}

#' CCplot
#' @description Plots the cell composition differences between groups.
#' @param tidyCC The list of tibbles returned by \code{DMRichR::CCstats()}
#' @param cellComposition Character vector of the cell composition dataset
#' @param testCovariate The factor to test for differences between groups
#' @param adjustCovariate The covariate(s) to adjust for between groups
#' @param matchCovariate Another covariate to adjust for between groups (for dmrseq compatibility)
#' @return A \code{ggplot} object that can be viewed by calling it, saved with \code{ggplot2::ggsave()},
#'  or further modified by adding \code{ggplot2} syntax.
#' @import tidyverse
#' @import ggsci
#' @export CCplot
CCplot <- function(tidyCC = tidyCC,
                   testCovariate = testCovariate,
                   adjustCovariate = NULL,
                   matchCovariate = NULL){
  
  tidyCC$summary %>%
    dplyr::select(one_of(!!testCovariate, !!adjustCovariate, !!matchCovariate, !!IDs)) %>% 
    dplyr::group_by_(testCovariate) %>%
    dplyr::summarise_at(IDs, mean) %>%
    tidyr::pivot_longer(cols = all_of(IDs),
                        names_to = "Cell Type",
                        values_to = "Cell Composition") %>% 
    ggplot(aes(fill = `Cell Type`,
               y = `Cell Composition`,
               x = !!rlang::sym(testCovariate))
    ) + 
    geom_bar(position = "stack",
             stat = "identity") +
    scale_y_continuous(expand = c(0, 0)) +
    ggsci::scale_fill_aaas() +
    scale_x_discrete(labels = function(x){stringr::str_wrap(x, width = 10)}) +
    theme_classic(base_size = 14)
}
  