#' bsseqLift
#' @description LiftOver a hg38 bsseq objet to hg19 coordinates
#' @param bs.filtered.bsseq A \code{bsseq} object with hg38 coordinates
#' @return A \code{bsseq} object with hg19 coordinates
#' @importFrom magrittr %>%
#' @importFrom rtracklayer liftOver
#' @importFrom AnnotationHub AnnotationHub
#' @import GenomicRanges
#' @importFrom glue glue
#' @importClassesFrom bsseq BSseq 
#' @export bsseqLift
#' 
bsseqLift <- function(bs.filtered.bsseq = bs.filtered.bsseq){
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

#' prepareCC
#' @description Prepare a \code{bsseq} object for cell composition estimation
#' @param bs.filtered.bsseq A \code{bsseq} object
#' @param genome Character string of genome symbol c("hg38", "hg19").
#' @return A \code{bsseq} object with hg19 coordinates
#' @importFrom magrittr %>%
#' @importFrom GenomeInfoDb dropSeqlevels
#' @importClassesFrom bsseq BSseq 
#' @export prepareCC
#' 
prepareCC <- function(bs.filtered.bsseq = bs.filtered.bsseq,
                      genome = genome){
  
  stopifnot(genome %in% c("hg38", "hg19"))
  
  if(genome == "hg38"){
    bs.filtered.bsseq.cc <- bsseqLift(bs.filtered.bsseq)
  }else{
    bs.filtered.bsseq.cc <- bs.filtered.bsseq
  }

  # Keep only autosomes
  bs.filtered.bsseq.cc <- bs.filtered.bsseq.cc %>%
    GenomeInfoDb::dropSeqlevels(c("chrX", "chrY", "chrM"), pruning.mode = "coarse") %>%
    unique()
}

#' arrayRanges
#' @description Obtain hg19 EPIC array coordinates
#' @return A \code{GRanges} object of hg19 EPIC coordinates
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr select
#' @importFrom GenomicRanges makeGRangesFromDataFrame granges
#' @importFrom magrittr %>%
#' @importFrom minfi getAnnotation
#' @references \url{https://support.bioconductor.org/p/78652/}
#' @export arrayRanges
#' 
arrayRanges <- function(){
  
  if(!require(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)){
    BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
    library("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
    }
  
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
    GenomicRanges::granges()
  
  return(array)
}

#' Houseman
#' @description Utilize the Houseman method to estimate whole-blood cell composition
#' @param bs.filtered.bsseq.cc A \code{bsseq} object that has been lifted over to hg19
#' @return Cell composition estimates
#' @import GenomicRanges
#' @importFrom IRanges subsetByOverlaps
#' @importFrom magrittr %>%
#' @importFrom bsseq getMeth
#' @importClassesFrom bsseq BSseq 
#' @export Houseman
#' 
Houseman <- function(bs.filtered.bsseq.cc = bs.filtered.bsseq.cc){
  
  # Subset bsseq for sites on EPIC
  EPIC <- DMRichR::arrayRanges()
  
  bsseq.filtered.EPIC <- bs.filtered.bsseq.cc %>% 
    IRanges::subsetByOverlaps(EPIC)
  
  names(bsseq.filtered.EPIC) <- EPIC %>%
    IRanges::subsetByOverlaps(bs.filtered.bsseq.cc) %>% 
    names()
  
  # liftOver to EPIC
  testingBeta <- bsseq.filtered.EPIC %>%
    bsseq::getMeth()
  
  rownames(testingBeta) <- rownames(bsseq.filtered.EPIC)
  
  # Estimates
  if(!require(FlowSorted.Blood.EPIC)){
    BiocManager::install("Immunomethylomics/FlowSorted.Blood.EPIC")
    library(FlowSorted.Blood.EPIC)
  }
  
  testingBeta <- testingBeta[rownames(testingBeta) %in% FlowSorted.Blood.EPIC::IDOLOptimizedCpGs,]
  
  CC <- FlowSorted.Blood.EPIC::projectCellType_CP(testingBeta,
                                                  FlowSorted.Blood.EPIC::IDOLOptimizedCpGs.compTable)
  
  return(CC)
}

#' CCstats
#' @description Computes the cell composition differences between groups while adjusting for the provided covariates. 
#'  The differences are tested for using an ANOVA through the \code{\link[stats]{aov}} function.
#' @param CC A \code{matrix} from \code{DMRichR::Houseman()}
#' @param bs.filtered.bsseq.cc Smoothed \code{bsseq} object with design matrix in pData
#' @param testCovariate The factor to test for differences between groups
#' @param adjustCovariate The covariate(s) to adjust for between groups
#' @param matchCovariate Another covariate to adjust for between groups (for dmrseq compatibility)
#' @return A list of tibbles with the statsitics and the values used for the tests
#' @references \url{https://cran.r-project.org/web/packages/broom/vignettes/broom_and_dplyr.html}
#' @import lsmeans
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr as_tibble full_join select group_by summarise_at mutate
#' @importFrom GenomicRanges as.data.frame
#' @importFrom magrittr %>%
#' @importFrom tidyr pivot_longer nest
#' @importFrom purrr map
#' @importFrom broom tidy
#' @importClassesFrom bsseq BSseq
#' @importMethodsFrom bsseq pData
#' @export CCstats
#' 
CCstats <- function(CC = NULL,
                    bs.filtered.bsseq.cc = bs.filtered.bsseq.cc,
                    testCovariate = testCovariate,
                    adjustCovariate = NULL,
                    matchCovariate = NULL){
  
  # Tidy --------------------------------------------------------------------
  
  IDs <- c("Neu", "NK", "Bcell" , "CD4T", "CD8T", "Mono")
  
  tidyCC <- CC %>% 
    tibble::rownames_to_column("Sample") %>% 
    dplyr::as_tibble() %>% 
    dplyr::full_join(bs.filtered.bsseq.cc %>%
                       pData() %>%
                       as.data.frame() %>% 
                       tibble::rownames_to_column("Sample") %>%
                       dplyr::as_tibble(),
                     .)
  
  summary <- tidyCC %>%
    dplyr::select(one_of(!!testCovariate, !!adjustCovariate, !!matchCovariate, !!IDs)) %>% 
    dplyr::group_by(!!as.name(testCovariate)) %>%
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
  
  pairWise <- tidyCC %>%
    tidyr::pivot_longer(cols = all_of(IDs),
                        names_to = "cellType",
                        values_to = "cellCount") %>% 
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
#' @import ggplot2
#' @importFrom dplyr select group_by summarise_at
#' @importFrom tidyr pivot_longer
#' @importFrom ggsci scale_fill_aaas
#' @importFrom stringr str_wrap
#' @export CCplot
#' 
CCplot <- function(tidyCC = tidyCC,
                   testCovariate = testCovariate,
                   adjustCovariate = NULL,
                   matchCovariate = NULL){
  
  IDs <- c("Neu", "NK", "Bcell" , "CD4T", "CD8T", "Mono")
  
  tidyCC$summary %>%
    dplyr::select(one_of(!!testCovariate, !!adjustCovariate, !!matchCovariate, !!IDs)) %>% 
    dplyr::group_by(!!as.name(testCovariate)) %>%
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

#' @title Finding cell type specific differentially methylated regions
#' 
#' @description This function is modified from \code{methylCC}
#' to use flow sorted array reference methylomes for either blood, 
#' cord blood, or brain to identify cell type specific 
#' differentially methylated regions.  
#'
#' @param verbose TRUE/FALSE argument specifying if verbose
#' messages should be returned or not. Default is TRUE.
#' @param gr_target Default is NULL. However, the user 
#' can provide a GRanges object from the \code{object} 
#' in \code{estimatecc}. Before starting the procedure to 
#' find differentially methylated regions, the intersection
#' of the \code{gr_target} and GRanges object from the 
#' reference methylomes (\code{FlowSorted.Blood.450k}). 
#' @param include_cpgs TRUE/FALSE. Should individual CpGs
#' be returned. Default is FALSE. 
#' @param include_dmrs TRUE/FALSE. Should differentially 
#' methylated regions be returned. Default is TRUE. User
#' can turn this to FALSE and search for only CpGs. 
#' @param num_cpgs The max number of CpGs to return 
#' for each cell type. Default is 50.
#' @param num_regions The max number of DMRs to return 
#' for each cell type. Default is 50. 
#' @param bumphunter_beta_cutoff The \code{cutoff} threshold 
#' in \code{bumphunter()} in the \code{bumphunter} package. 
#' @param dmr_up_cutoff A cutoff threshold for identifying 
#' DMRs that are methylated in one cell type, but not in the 
#' other cell types. 
#' @param dmr_down_cutoff A cutoff threshold for identifying 
#' DMRs that are not methylated in one cell type, but 
#' methylated in the other cell types.
#' @param dmr_pval_cutoff  A cutoff threshold for the p-values 
#' when identifying DMRs that are methylated in one cell 
#' type, but not in the other cell types (or vice versa). 
#' @param cpg_pval_cutoff A cutoff threshold for the p-values 
#' when identifying differentially methylated CpGs that are
#' methylated in one cell type, but not in the other cell
#' types (or vice versa). 
#' @param cpg_up_dm_cutoff A cutoff threshold for identifying 
#' differentially methylated CpGs that are methylated in 
#' one cell type, but not in the other cell types. 
#' @param cpg_down_dm_cutoff A cutoff threshold for identifying 
#' differentially methylated CpGs that are not methylated in 
#' one cell type, but are methylated in the other cell types.  
#' @param pairwise_comparison TRUE/FAlSE of whether all pairwise
#' comparisons (e.g. methylated in Granulocytes and Monocytes, 
#' but not methylated in other cell types). Default if FALSE. 
#' @param mset_train_flow_sort Default is NULL and will select 
#' select \code{"FlowSorted.Blood.450k"}. Alternatively, the user
#' can specify another existing dataset \code{c("FlowSorted.Blood.EPIC",
#' "FlowSorted.CordBloodCombined.450k", "FlowSorted.CordTissueAndBlood.EPIC",
#' "FlowSorted.CordBlood.450k", "FlowSorted.DLPFC.450k")}.
#'  Additionally, a user can provide a 
#' \code{MethylSet} object after processing the 
#' \code{FlowSorted.Blood.450k} dataset. The default
#' normalization is \code{preprocessQuantile()}.
#' 
#' @return A list of data frames and GRanges objects.
#' 
#' @import GenomicRanges
#' @importFrom Biobase pData
#' @importFrom bumphunter clusterMaker loessByCluster bumphunter 
#' @importFrom genefilter rowttests
#' @importFrom plyranges arrange 
#' @importFrom S4Vectors queryHits
#' @importFrom stats model.matrix
#' @importFrom utils head
#' @importFrom glue glue
#' @importFrom magrittr %>%
#' @importFrom minfi preprocessIllumina preprocessQuantile preprocessNoob mapToGenome
#' @export find_dmrs2
#' 
find_dmrs2 <- function(verbose = TRUE, gr_target = NULL,
                       include_cpgs = FALSE, include_dmrs = TRUE,
                       num_cpgs = 50, num_regions = 50, 
                       bumphunter_beta_cutoff = 0.2, 
                       dmr_up_cutoff = 0.5, dmr_down_cutoff = 0.4,
                       dmr_pval_cutoff = 1e-11, cpg_pval_cutoff = 1e-08,
                       cpg_up_dm_cutoff = 0, cpg_down_dm_cutoff = 0, 
                       pairwise_comparison = FALSE,
                       mset_train_flow_sort = NULL) {
  
  if(!require(methylCC)){
    BiocManager::install("methylCC")
    library("methylCC")
    }
  
  print(glue::glue("Finding hg19 cell type specific DMRs using {mset_train_flow_sort}"))
  
  if(mset_train_flow_sort == "FlowSorted.Blood.450k"){
    
    if(!require(FlowSorted.Blood.450k)){
      BiocManager::install("FlowSorted.Blood.450k")
      library("FlowSorted.Blood.450k")
      }
    data(FlowSorted.Blood.450k)
    
    dataset <- mset_train_flow_sort
    mset_train_flow_sort <- FlowSorted.Blood.450k
    mset_train_flow_sort <- mset_train_flow_sort %>%
      updateObject() %>%
      preprocessQuantile() %>% 
      mapToGenome(mergeManifest = FALSE)
    
    IDs <- c("Gran", "CD4T", "CD8T", "Bcell", "Mono", "NK")
    
    # remove outliers
    mset_train_flow_sort <- mset_train_flow_sort[, 
                                                 pData(mset_train_flow_sort)$Sample_Name != "CD8+_105"]
    
  }else if(mset_train_flow_sort == "FlowSorted.Blood.EPIC"){
    
    if(!require(FlowSorted.Blood.EPIC)){
      BiocManager::install("FlowSorted.Blood.EPIC")}
    if(!require(ExperimentHub)){
      BiocManager::install("ExperimentHub")}
    
    dataset <- mset_train_flow_sort
    mset_train_flow_sort <- ExperimentHub::ExperimentHub()[["EH1136"]] %>%
      updateObject() %>%
      preprocessQuantile() %>%
      mapToGenome(mergeManifest = FALSE)
    
    IDs <- c("Neu", "NK", "Bcell" , "CD4T", "CD8T", "Mono")
    
  }else if(mset_train_flow_sort == "FlowSorted.CordTissueAndBlood.EPIC"){
    
    if(!require(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)){
      BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b2.hg19")}
    if(!require(FlowSorted.CordTissueAndBlood.EPIC)){
      install.packages(paste0("https://karnanilab.com/Tools/",
                              "FlowSorted.CordTissueAndBlood.EPIC/",
                              "FlowSorted.CordTissueAndBlood.EPIC_1.0.1.tar.gz"),
                       repos = NULL, type = "source")
      library(FlowSorted.CordTissueAndBlood.EPIC)
      }
    data(FlowSorted.CordTissueAndBlood.EPIC)
    
    dataset <- mset_train_flow_sort
    mset_train_flow_sort <- FlowSorted.CordTissueAndBlood.EPIC %>% 
      updateObject() %>%
      preprocessNoob() %>%
      mapToGenome(mergeManifest = FALSE)
    
    IDs <- c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran")
    
  }else if(mset_train_flow_sort == "FlowSorted.CordBloodCombined.450k"){
    
    if(!require(FlowSorted.CordBloodCombined.450k)){
      BiocManager::install("FlowSorted.CordBloodCombined.450k")}
    if(!require(ExperimentHub)){
      BiocManager::install("ExperimentHub")}
    
    dataset <- mset_train_flow_sort
    mset_train_flow_sort <- ExperimentHub::ExperimentHub()[["EH2256"]] %>%
      updateObject() %>%
      preprocessNoob() %>%
      mapToGenome(mergeManifest = FALSE)
    
    IDs <- c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran", "nRBC")
  
  }else if(database == "FlowSorted.CordBlood.450k"){
    if(!require(FlowSorted.CordBlood.450k)){
      BiocManager::install("FlowSorted.CordBlood.450k")
      library("FlowSorted.CordBlood.450k")
      }
    data(FlowSorted.CordBlood.450k)
    
    dataset <- mset_train_flow_sort
    mset_train_flow_sort <- FlowSorted.CordBlood.450k
    mset_train_flow_sort <- mset_train_flow_sort %>%
      updateObject() %>%
      preprocessNoob() %>%
      mapToGenome(mergeManifest = FALSE)
    
  }else if(mset_train_flow_sort == "FlowSorted.DLPFC.450k"){
    
    if(!require(FlowSorted.DLPFC.450k)){
      BiocManager::install("FlowSorted.DLPFC.450k")
      library("FlowSorted.DLPFC.450k")
      }
    
    dataset <- mset_train_flow_sort
    mset_train_flow_sort <- FlowSorted.DLPFC.450k
    mset_train_flow_sort <- mset_train_flow_sort %>%
      updateObject() %>%
      preprocessQuantile() %>%
      mapToGenome(mergeManifest = FALSE)
    
    IDs <- c("NeuN_pos", "NeuN_neg")
    
  }else{
    # remove outliers
    dataset <- "FlowSorted.Blood.450k"
    mset_train_flow_sort <- mset_train_flow_sort[, 
                                                 pData(mset_train_flow_sort)$Sample_Name != "CD8+_105"]
    IDs <- c("Gran", "CD4T", "CD8T", "Bcell", "Mono", "NK")
  }
  
  # create training object to identify DMRs
  mset_train_flow_sort <- mset_train_flow_sort[, 
                                               (pData(mset_train_flow_sort)$CellType %in% IDs) ]
  
  # find celltype specific regions using only overlap CpGs in target object
  if(!is.null(gr_target)){ 
    # which of the 450K CpGs overlap with the target CpGs
    zz <- findOverlaps(granges(mset_train_flow_sort), gr_target) 
    mset_train_flow_sort <- mset_train_flow_sort[queryHits(zz), ]
    if(verbose){
      mes <- "[estimatecc] gr_target is not null. Using %s overlapping CpGs."
      message(sprintf(mes, nrow(mset_train_flow_sort)))
    }
  }  
  
  # extract beta values, phenotypic information and GRanges objects
  pd <- as.data.frame(pData(mset_train_flow_sort))
  gr <- granges(mset_train_flow_sort)
  p_beta <- getBeta(mset_train_flow_sort) # beta values # , type = "Illumina"
  
  if(dataset == "FlowSorted.Blood.450k"){
    colnames(p_beta) = pd$Sample_Name = rownames(pd) = 
      gsub("\\+","", pd$Sample_Name)
  }else if(dataset == "FlowSorted.Blood.EPIC" |
           dataset == "FlowSorted.CordBlood.EPIC" |
           dataset == "FlowSorted.CordBloodCombined.450k" |
           dataset == "FlowSorted.CordBlood.450k"){
    colnames(p_beta) = pd$Sample_Name = rownames(pd) =
      paste(pd$CellType, pd$Sample_Name, pd$Slide, sep="_")
  }else if(dataset == "FlowSorted.DLPFC.450k"){
    colnames(p_beta) = pd$Sample_Name = rownames(pd) =
      paste(pd$CellType, pd$SampleID, sep="_")
  }
  
  cell <- factor(pd$CellType, levels = IDs)
  cell_levels <- levels(cell)
  
  # extract chromosome and position information for each probe in 
  #   450k array (need this for regions)
  chr <- as.character(seqnames(gr))
  pos <- start(gr)
  cl <- clusterMaker(chr, pos) # Create clusters using clusterMaker()
  
  # define design matrix to search for DMRs
  xmat = cbind(rep(1, length(cell)), model.matrix(~cell - 1))
  colnames(xmat) = c("Intercept", cell_levels)
  
  if(pairwise_comparison){
    all_poss = as.matrix(expand.grid(c(0,1), c(0,1), c(0,1), 
                                     c(0,1), c(0,1), c(0,1)))
    # remove the cases containing all methylated or unmethylated. 
    all_poss = all_poss[2:32,] 
    all_poss <- (all_poss == TRUE)
    colnames(all_poss) <- cell_levels
  } else { 
    all_poss = diag(length(cell_levels))
    all_poss <- (all_poss == TRUE)
    colnames(all_poss) <- cell_levels
  }
  
  regions_all <- GRanges() 
  zmat <- c() # regions_all, will contain all celltype-specific DMRs
  for(ind in seq_len(nrow(all_poss))){
    if(verbose){
      if(include_dmrs & include_cpgs){
        mes <- "[estimatecc] Searching for %s cell type-specific 
                regions and CpGs."
        message(sprintf(mes, paste(cell_levels[all_poss[ind,]], collapse=",")))
      } 
      if(include_dmrs & !include_cpgs) {
        mes <- "[estimatecc] Searching for %s cell type-specific regions."
        message(sprintf(mes, paste(cell_levels[all_poss[ind,]], collapse=",")))
      }
      if(!include_dmrs & include_cpgs) {
        mes <- "[estimatecc] Searching for %s cell type-specific CpGs"
        message(sprintf(mes, paste(cell_levels[all_poss[ind,]], collapse=",")))
      }
    }
    
    x_ind = cbind("Intercept" = xmat[, "Intercept"],
                  "cellTypes" = rowSums(
                    as.matrix(xmat[, cell_levels[all_poss[ind,]] ],
                              ncols = length(cell_levels[all_poss[ind,]]))))
    
    if(!include_dmrs){ 
      gr_regions_up <- GRanges()
      gr_regions_down <- GRanges()
    }
    
    if(include_dmrs){
      bumps = bumphunter(object = p_beta, design = x_ind, 
                         chr = chr, pos = pos, cluster = cl, 
                         cutoff = bumphunter_beta_cutoff, 
                         B = 0, smooth = FALSE, 
                         smoothFunction = loessByCluster)
      
      # y_regions are the beta values collapsed (CpGs averaged) by regions 
      # from bumphunter
      y_regions <- t(apply(bumps$table[,7:8], 1, function(z){
        colMeans(p_beta[(z[1]):(z[2]),,drop=FALSE]) } ))
      
      tmp <- genefilter::rowttests(y_regions,factor(x_ind[,"cellTypes"]))
      bumps$table$p.value <- tmp$p.value
      bumps$table$dm <- tmp$dm 
      bumps$table$dmr_up_max_diff <- 
        apply(abs(sweep(y_regions, 2, x_ind[,"cellTypes"], FUN = "-")), 1, max)
      bumps$table$dmr_down_max_diff <- 
        apply(abs(sweep(y_regions, 2, (1 - x_ind[,"cellTypes"]), FUN = "-")), 
              1, max)
      
      # # Only include region with more than 1 CpG (L > 1)
      # #       OR only 1 CpG in region if no other larger regions possible
      L = dm <- NULL 
      keep_ind_regions <- (bumps$table$L > 1 | 
                             (bumps$table$L==1 & bumps$table$clusterL == 1)) & 
        (bumps$table$p.value < dmr_pval_cutoff)  # ideally less than 1e-11
      
      bump_mat_up <- bumps$table[keep_ind_regions & bumps$table$dm < 0 &
                                   # ideally less than 0.6
                                   bumps$table$dmr_up_max_diff<dmr_up_cutoff,] 
      bump_mat_up <- bump_mat_up[order(-bump_mat_up$L, bump_mat_up$dm), ]
      if(nrow(bump_mat_up) > 0){
        gr_regions_up <- makeGRangesFromDataFrame(bump_mat_up, 
                                                  keep.extra.columns=TRUE)
        mcols(gr_regions_up)$dmr_status <- rep("DMR", length(gr_regions_up))
        gr_regions_up <- gr_regions_up[, names(mcols(gr_regions_up)) %in% 
                                         c("indexStart", "indexEnd", "L", "dm", "p.value", 
                                           "dmr_status", "dmr_up_max_diff")]
        names(mcols(gr_regions_up))[
          names(mcols(gr_regions_up)) == "dmr_up_max_diff"] <- "dmr_max_diff"
        
        gr_regions_up <- gr_regions_up %>% arrange(-L, dm) %>% 
          head(num_regions)
      } else {
        gr_regions_up <- GRanges()
      }
      
      bump_mat_down <- bumps$table[(keep_ind_regions) & bumps$table$dm > 0 & 
                                     bumps$table$dmr_down_max_diff < 
                                     dmr_down_cutoff,] # ideally less than 0.8
      bump_mat_down <- bump_mat_down[order(-bump_mat_down$L, 
                                           -bump_mat_down$dm), ]
      if(nrow(bump_mat_down) > 0){
        gr_regions_down <- makeGRangesFromDataFrame(bump_mat_down, 
                                                    keep.extra.columns=TRUE)
        mcols(gr_regions_down)$dmr_status <- 
          rep("DMR", length(gr_regions_down))
        gr_regions_down <- gr_regions_down[, names(mcols(gr_regions_down)) %in%
                                             c("indexStart", "indexEnd", "L", "dm", "p.value", 
                                               "dmr_status", "dmr_down_max_diff")]
        names(mcols(gr_regions_down))[names(mcols(gr_regions_down)) 
                                      == "dmr_down_max_diff"] <- "dmr_max_diff"
        gr_regions_down <- gr_regions_down %>% 
          arrange(-L, -dm) %>% 
          head(num_regions)
      } else {
        gr_regions_down <- GRanges()
      }
    }  
    
    if(include_cpgs){
      tstats <- genefilter::rowttests(p_beta, factor(x_ind[,"cellTypes"]))
      tstats <- tstats[(tstats[, "p.value"] < cpg_pval_cutoff),] 
      
      tstats_up <- tstats[order(tstats[, "dm"], decreasing = FALSE), ]
      # at a min should be less than 0
      tstats_up <- tstats_up[tstats_up$dm < cpg_up_dm_cutoff,] 
      
      probe_keep <- rownames(tstats_up)[seq_len(min(nrow(tstats_up), 
                                                    num_cpgs))]
      if(length(probe_keep) > 0){
        gr_probe <- granges(mset_train_flow_sort[probe_keep,])
        mcols(gr_probe) <- tstats[probe_keep, c("dm", "p.value")]
        mcols(gr_probe)$L <- rep(1, length(probe_keep))
        mcols(gr_probe)$indexStart <- match(probe_keep,
                                            rownames(mset_train_flow_sort))
        mcols(gr_probe)$indexEnd <- match(probe_keep, 
                                          rownames(mset_train_flow_sort))
        mcols(gr_probe)$dmr_status <- rep("CpG", length(gr_probe))
        gr_regions_up <- unique(c(gr_regions_up, 
                                  gr_probe[,c("indexStart", "indexEnd", 
                                              "L", "p.value", "dm", 
                                              "dmr_status")]))
        gr_regions_up <- gr_regions_up %>% arrange(-L, dm) %>% 
          head(num_regions)
      } 
      
      tstats_down <- tstats[order(tstats[, "dm"], decreasing = TRUE), ]
      # at a min should be greater than 0
      tstats_down <- tstats_down[tstats_down$dm > cpg_down_dm_cutoff,] 
      probe_keep <- rownames(tstats_down)[seq_len(min(nrow(tstats_down), 
                                                      num_cpgs))]
      if(length(probe_keep) > 0){
        gr_probe <- granges(mset_train_flow_sort[probe_keep,])
        mcols(gr_probe) <- tstats[probe_keep, c("dm", "p.value")]
        mcols(gr_probe)$L <- rep(1, length(probe_keep))
        mcols(gr_probe)$indexStart <- match(probe_keep,
                                            rownames(mset_train_flow_sort))
        mcols(gr_probe)$indexEnd <- match(probe_keep, 
                                          rownames(mset_train_flow_sort))
        mcols(gr_probe)$dmr_status <- rep("CpG", length(gr_probe))
        gr_regions_down <- unique(c(gr_regions_down, 
                                    gr_probe[,c("indexStart", "indexEnd", 
                                                "L", "p.value", "dm", 
                                                "dmr_status")]))
        gr_regions_down <- gr_regions_down %>% 
          arrange(-L, -dm) %>% 
          head(num_regions)
      } 
    }
    
    mcols(gr_regions_up)$status <- rep("Up", length(gr_regions_up))
    mcols(gr_regions_down)$status <- rep("Down", length(gr_regions_down))
    bump_mat_all <- c(gr_regions_up, gr_regions_down)
    mcols(bump_mat_all)$cellType <- rep(paste(cell_levels[all_poss[ind,]], 
                                              collapse=","), 
                                        length(bump_mat_all)) 
    
    if(verbose){
      if(include_dmrs & include_cpgs){
        mes <- "[estimatecc] Found %s %s cell type-specific regions and CpGs."
        message(sprintf(mes, length(bump_mat_all), 
                        paste(cell_levels[all_poss[ind,]], collapse=",")))
      } 
      if(include_dmrs & !include_cpgs) {
        mes <- "[estimatecc] Found %s %s cell type-specific regions."
        message(sprintf(mes, length(bump_mat_all), 
                        paste(cell_levels[all_poss[ind,]], collapse=",")))
      }
      if(!include_dmrs & include_cpgs) {
        mes <- "[estimatecc] Found %s %s cell type-specific CpGs."
        message(sprintf(mes, length(bump_mat_all),
                        paste(cell_levels[all_poss[ind,]], collapse=",")))
      }
    }
    
    if(length(bump_mat_all) > 0){ 
      regions_all <- c(regions_all, bump_mat_all)
    }
    if(length(gr_regions_up) > 0){
      zmat <- rbind(zmat, t(replicate(min(length(gr_regions_up), num_regions), 
                                      as.numeric(all_poss[ind,]))))
    }
    if(length(gr_regions_down) > 0){
      zmat <- rbind(zmat, t(replicate(min(length(gr_regions_down), 
                                          num_regions), 
                                      as.numeric(!all_poss[ind,]))))
    }
  }
  colnames(zmat) <- cell_levels
  
  y_regions <- t(apply(
    as.data.frame(mcols(regions_all))[,seq_len(2)],1,function(ind){
      colMeans(p_beta[(ind[1]):(ind[2]),,drop=FALSE])
    }))
  
  profiles <- vapply(methylCC:::.splitit(cell), 
                     FUN = function(ind){ rowMeans(y_regions[,ind])}, 
                     FUN.VALUE = numeric(nrow(y_regions)))
  
  removeMe <- duplicated(regions_all)
  list(regions_all = regions_all[!removeMe,], 
       zmat = zmat[!removeMe,], 
       y_regions = y_regions[!removeMe,], 
       profiles = profiles[!removeMe,],
       cell = cell, cell_mat = all_poss, 
       cell_levels = cell_levels, pd = pd)
} 
  