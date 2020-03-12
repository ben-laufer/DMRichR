#' processBismark
#' @description Process bismark cytosine reports into bsseq objects with design matrix pData
#' @param files List of cytosine report file paths
#' @param meta Design matrix table with sample name in the Name column 
#' @param testCovar Factor of interest (testCovariate)
#' @param adjustCovar Variables to adjust for (adjustCovariate)
#' @param matchCovar Variable to block for when constructing permutations (matchCovariate)
#' @param Cov CpG coverage cutoff (1x recommended)
#' @param mc.cores Number of cores to use
#' @param per.Group Percent of samples per a group to apply the CpG coverage cutoff to
#' @import bsseq
#' @importFrom openxlsx read.xlsx
#' @import tidyverse
#' @import optparse
#' @importFrom parallel mclapply
#' @importFrom glue glue
#' @export processBismark
processBismark <- function(files = list.files(path = getwd(), pattern = "*.txt.gz"),
                           meta = openxlsx::read.xlsx("sample_info.xlsx", colNames = TRUE) %>% dplyr::mutate_if(is.character, as.factor),
                           testCovar = testCovariate,
                           adjustCovar = NULL,
                           matchCovar = NULL,
                           Cov = coverage,
                           mc.cores = cores,
                           per.Group = perGroup){
  
  cat("\n[DMRichR] Processing Bismark cytosine reports \t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
  start_time <- Sys.time()
  print(glue::glue("Selecting files..."))
  files.idx <- pmatch(meta$Name, files)
  files <- files[files.idx]
  #names <- as.data.frame(gsub( "_.*$","", files[files.idx])) # For colData, but jumbles file order with parallel processing
  #colnames(names) <- "Name"
  #rownames(names) <- names[,1]
  #names[,1] <- NULL
  
  # glue::glue("Determining parallelization...") # Does not work on some clusters due to use of BiocParallel, but speeds up desktops 
  # if(mc.cores >= 4){
  #  BPPARAM <- BiocParallel::MulticoreParam(workers = floor(mc.cores/4), progressbar = TRUE)
  #  nThread <- as.integer(floor(mc.cores/floor(mc.cores/4)))
  #  glue::glue("Parallel processing will be used with {floor(mc.cores/4)} cores consisting of {nThread} threads each")
  # }else if(mc.cores < 4){
  #  BPPARAM <- BiocParallel::MulticoreParam(workers = 1, progressbar = TRUE)
  #  nThread <- as.integer(1)
  #  glue::glue("Parallel processing will not be used")
  # }
  
  print(glue::glue("Reading cytosine reports..."))
  bs <- read.bismark(files = files,
                     #colData = names,
                     rmZeroCov = FALSE,
                     strandCollapse = TRUE,
                     verbose = TRUE,
                     BPPARAM = MulticoreParam(workers = mc.cores, progressbar = FALSE), # BPPARAM # bpparam() # MulticoreParam(workers = mc.cores, progressbar = TRUE)
                     nThread = 1) # 1L # nThread
  
  print(glue::glue("Assigning sample metadata with {testCovar} as factor of interest..."))
  sampleNames(bs) <- gsub( "_.*$","", sampleNames(bs))
  meta <- meta[order(match(meta[,1],sampleNames(bs))),]
  stopifnot(sampleNames(bs) == as.character(meta$Name))
  pData(bs) <- cbind(pData(bs), meta[2:length(meta)])
  print(pData(bs))
  
  print(glue::glue("Filtering CpGs for {testCovar}..."))
  bs <- GenomeInfoDb::keepStandardChromosomes(bs, pruning.mode = "coarse")
  pData(bs)[[testCovar]] <- as.factor(pData(bs)[[testCovar]])
  loci.cov <- getCoverage(bs, type = "Cov")
  
  if(!is.null(adjustCovar)){
    excludeCovar <- NULL
    for(i in 1:length(adjustCovar)){
      if(is.numeric(pData(bs)[, adjustCovar[i]])){
        print(glue::glue("Assuming adjustment covariate {adjustCovar[i]} is continuous and excluding it from filtering..."))
        excludeCovar <- c(excludeCovar, adjustCovar[i])
        
      }else{
        print(glue::glue("Assuming adjustment covariate {adjustCovar[i]} is discrete and including it for filtering..."))
      }
    }
    adjustCovar <- adjustCovar[!adjustCovar %in% excludeCovar]
  }
  
  if(!is.null(matchCovar)){
    if(length(matchCovar) > 1){
      stop(print(glue::glue("Only one matching covariate can be used")))
      
    }else if(is.numeric(pData(bs)[, matchCovar])){
      stop(print(glue::glue("Matching covariate {matchCovar} must be discrete")))
      
    }else{
      print(glue::glue("Assuming matching covariate {matchCovar} is discrete and including it for filtering..."))
    }
  }
  
  covar.groups <- apply(pData(bs)[, as.character(c(testCovar, adjustCovar, matchCovar))] %>% as.data.frame(), 
                        MARGIN = 1, FUN = paste, collapse = "_") %>% 
    as.factor() # Covariate combination groups
  
  group.samples <- split(t(loci.cov >= Cov) %>% as.data.frame(), f = covar.groups) %>% 
    parallel::mclapply(FUN = as.matrix, mc.cores = mc.cores) %>% 
    parallel::mclapply(FUN = DelayedMatrixStats::colSums2, mc.cores = mc.cores) %>%
    simplify2array() %>%
    as.data.frame() # Samples in each cov.group meeting coverage threshold by CpG (slow)
  
  print(glue::glue("Making coverage filter table..."))
  per.Group.seq <- seq(0,1,0.05)
  covFilter <- NULL
  for(i in 1:length(per.Group.seq)){
    groups.n <- (table(covar.groups) * per.Group.seq[i]) %>% ceiling() %>% as.integer()
    per.Group.seq.test <- mapply(function(x, y){x >= y}, 
                                 x = group.samples, 
                                 y = (table(covar.groups) * per.Group.seq[i]) %>%
                                   ceiling() %>%
                                   as.integer()) # Test if enough samples are in each group by CpG
    CpGs <- sum(DelayedMatrixStats::rowSums2(per.Group.seq.test) >= length(unique(covar.groups))) # Total CpGs meeting coverage threshold in at least per.Group of all covariate combos
    temp <- c(per.Group.seq[i] * 100, groups.n, CpGs, round(CpGs * 100 / length(bs), 2))
    covFilter <- rbind(covFilter, temp)
  }
  covFilter <- as.data.frame(covFilter, row.names = 1:nrow(covFilter))
  colnames(covFilter) <- c("perGroup", paste("n", levels(covar.groups), sep = "_"), "nCpG", "perCpG")
  print(covFilter)
  
  if(per.Group == 1){
    print(glue::glue("Filtering for {Cov}x coverage in all samples"))
    sample.idx <- which(pData(bs)[[testCovar]] %in% levels(pData(bs)[[testCovar]]))
    loci.idx <- which(DelayedMatrixStats::rowSums2(getCoverage(bs, type = "Cov") >= Cov) >= length(sample.idx))
    bs.filtered <- bs[loci.idx, sample.idx]
  
  }else if(per.Group < 1){
    print(glue::glue("Filtering for {Cov}x coverage in at least {per.Group*100}% of samples for \\
                                 all combinations of covariates..."))
    sample.idx <- which(pData(bs)[[testCovar]] %in% levels(pData(bs)[[testCovar]]))
    per.Group.test <- mapply(function(x, y){x >= y}, 
                             x = group.samples, 
                             y = (table(covar.groups) * per.Group) %>% ceiling() %>% as.integer()) # Test if enough samples are in each group by CpG
    loci.idx <- which(DelayedMatrixStats::rowSums2(per.Group.test) >= length(unique(covar.groups))) # Which CpGs meet coverage threshold in at least per.Group of all covariate combos
    bs.filtered <- bs[loci.idx, sample.idx]
      
  }else if(per.Group > 1){
    stop(print(glue::glue("perGroup is {per.Group} and cannot be greater than 1, which is 100% of samples")))
    
  }else{
    stop(print(glue::glue("processBismark arguments")))
  } 
  
  print(glue::glue("processBismark timing..."))
  end_time <- Sys.time()
  print(end_time - start_time)
  
  print(glue::glue("Before filtering for {Cov}x coverage there were {nrow(bs)} CpGs, \\
                         after filtering there are {nrow(bs.filtered)} CpGs, \\
                         which is {round(nrow(bs.filtered)/nrow(bs)*100,1)}% of all CpGs."))
  
  return(bs.filtered)
}
