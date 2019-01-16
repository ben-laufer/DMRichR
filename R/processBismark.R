#' processBismark
#' @description Process bismark cytosine reports into bsseq objects with design matrix pData
#' @param files List of cytosine report file paths
#' @param meta Design matrix table with sample name in the Name column 
#' @param groups Factor of interest (testCovariate)
#' @param Cov CpG coverage cutoff (1x recommended)
#' @param mc.cores Number of cores to use
#' @param per.Group Percent of samples per a group to apply the CpG coverage cutoff to (only works for two factor testCovariates)
#' @import bsseq
#' @import openxlsx
#' @import tidyverse
#' @export processBismark
processBismark <- function(files = list.files(path = getwd(), pattern = "*.txt.gz"),
                           meta = read.xlsx("sample_info.xlsx", colNames = TRUE) %>% mutate_if(is.character,as.factor),
                           groups = testCovariate,
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
                     BPPARAM = MulticoreParam(workers = mc.cores, progressbar = TRUE), # BPPARAM # bpparam() # MulticoreParam(workers = mc.cores, progressbar = TRUE)
                     nThread = 1) # 1L # nThread
  
  print(glue::glue("Assigning sample metadata with {groups} as factor of interest..."))
  sampleNames(bs) <- gsub( "_.*$","", sampleNames(bs))
  meta <- meta[order(match(meta[,1],sampleNames(bs))),]
  stopifnot(sampleNames(bs) == as.character(meta$Name))
  pData(bs) <- cbind(pData(bs), meta[2:length(meta)])
  print(pData(bs))
  
  print(glue::glue("Filtering CpGs..."))
  bs <- GenomeInfoDb::keepStandardChromosomes(bs, pruning.mode = "coarse")
  pData(bs)[[groups]] <- as.factor(pData(bs)[[groups]])
  
  if(per.Group == 1){
    sample.idx <- which(pData(bs)[[groups]] %in% levels(pData(bs)[[groups]]))
    loci.idx <- which(DelayedMatrixStats::rowSums2(getCoverage(bs, type="Cov") >= Cov) >= length(sample.idx))
    bs.filtered <- bs[loci.idx, sample.idx]
  
  }else if(length(levels(pData(bs)[[groups]])) == 2 & per.Group < 1){
    print(glue::glue("Filtering for {Cov}x coverage in at least {per.Group*100}% of \\
                     {levels(pData(bs)[[groups]])[2]} and {levels(pData(bs)[[groups]])[1]} samples"))
    sample.idx <- which(pData(bs)[[groups]] %in% levels(pData(bs)[[groups]]))
    loci.cov <- getCoverage(bs, type = "Cov")
    ctrl.idx <- pData(bs)[[groups]] == levels(pData(bs)[[groups]])[1]
    exp.idx <- pData(bs)[[groups]] == levels(pData(bs)[[groups]])[2]
    loci.idx <- which(DelayedMatrixStats::rowSums2(loci.cov[, ctrl.idx] >= Cov) >= ceiling(per.Group * sum(ctrl.idx)) & 
                        DelayedMatrixStats::rowSums2(loci.cov[, exp.idx] >= Cov) >= ceiling(per.Group * sum(exp.idx))) 
    
    bs.filtered <- bs[loci.idx, sample.idx]
  
  }else if(per.Group > 1){
    stop(print(glue::glue("perGroup is {per.Group} and cannot be greater than 1, which is 100% of samples")))
    
  }else if(length(levels(pData(bs)[[groups]])) != 2 & per.Group < 1){
    stop(print(glue::glue("Samples cannot currently be filtered using a perGroup cutoff of {per.Group} \\
                          for a {length(levels(pData(bs)[[groups]]))} level testCovariate of {testCovariate}")))
  
  }else{
    stop(print(glue::glue("processBismark arguments")))
  } 
    
  print(glue::glue("processBismark timing..."))
  end_time <- Sys.time()
  print(end_time - start_time)
  
  print(glue::glue("Before filtering for {Cov}x coverage there were {nrow(bs)} CpGs, \\
                   after filtering for {Cov}x coverage there are {nrow(bs.filtered)} CpGs assayed"))
  
  return(bs.filtered)
}
