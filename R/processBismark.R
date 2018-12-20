#' processBismark
#' @description Process bismark cytosine reports into bsseq objects with design matrix pData
#' @param files List of cytosine report file paths
#' @param meta Design matrix table with sample name in the Name column 
#' @param groups Factor of interest (testCovariate)
#' @param Cov Coverage cutoff (1x recommended)
#' @param mc.cores Number of cores to use
#' @param per.ctrl Percent of control samples with coverage cutoff
#' @param per.exp Percent of experimental samples with coverage cutoff
#' @import bsseq
#' @import openxlsx
#' @import tidyverse
#' @export processBismark
processBismark <- function(files = list.files(path = getwd(), pattern = "*.txt.gz"),
                           meta = read.xlsx("sample_info.xlsx", colNames = TRUE) %>% mutate_if(is.character,as.factor),
                           groups = testCovariate,
                           Cov = coverage,
                           mc.cores = cores,
                           per.ctrl = 1,
                           per.exp = 1){
  cat("\n[DMRichR] Processing Bismark cytosine reports \t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
  start_time <- Sys.time()
  glue::glue("Selecting files...")
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
 
  glue::glue("Reading cytosine reports...")
  bs <- read.bismark(files = files,
                     #colData = names,
                     rmZeroCov = FALSE,
                     strandCollapse = TRUE,
                     verbose = TRUE,
                     BPPARAM = MulticoreParam(workers = mc.cores, progressbar = TRUE), # BPPARAM # bpparam() # MulticoreParam(workers = mc.cores, progressbar = TRUE)
                     nThread = 1) # 1L # nThread
  
  glue::glue("\n","Assigning sample metadata...")
  sampleNames(bs) <- gsub( "_.*$","", sampleNames(bs))
  meta <- meta[order(match(meta[,1],sampleNames(bs))),]
  stopifnot(sampleNames(bs) == as.character(meta$Name))
  pData(bs) <- cbind(pData(bs), meta[2:length(meta)])
  print(pData(bs))
  
  #glue::glue("\n", paste("Before filtering CpGs for ", Cov, "x coverage...", sep =""))
  glue::glue("\n","Filtering CpGs...")
  bs <- GenomeInfoDb::keepStandardChromosomes(bs, pruning.mode = "coarse")
  #print(head(getCoverage(bs, type = "Cov")))
  #print(bs)
  pData(bs)[[groups]] <- as.factor(pData(bs)[[groups]])
  saveRDS(bs, "unfiltered_BSseq_object.rds")
  
  stopifnot(length(levels(pData(bs)[[groups]])) == 2) # Filtering only works for two-group comparison
  stopifnot(per.ctrl <= 1 & per.exp <= 1 ) # per.ctrl and per.exp must be a decimal or 1
  sample.idx <- which(pData(bs)[[groups]] %in% levels(pData(bs)[[groups]]))
  #loci.idx <- which(DelayedMatrixStats::rowSums2(getCoverage(bs, type="Cov") >= Cov) >= length(sample.idx))
  #bs.filtered <- bs[loci.idx, sample.idx]
  sample.ctrl <- pData(bs)[[groups]] == levels(pData(bs)[[groups]])[1]
  sample.exp <- pData(bs)[[groups]] == levels(pData(bs)[[groups]])[2]
  loci.cov <- getCoverage(bs, type = "Cov") >= Cov
  loci.ctrl <- DelayedMatrixStats::rowSums2(loci.cov[,sample.ctrl]) >= ceiling(per.ctrl * sum(sample.ctrl))
  loci.exp <- DelayedMatrixStats::rowSums2(loci.cov[,sample.exp]) >= ceiling(per.exp * sum(sample.exp))
  bs.filtered <- bs[which(loci.ctrl & loci.exp), sample.idx]

  #glue::glue("\n", paste("After filtering CpGs for ", Cov, "x coverage...", sep =""))
  #print(head(getCoverage(bs.filtered, type = "Cov")))
  #print(bs.filtered)

  glue::glue("\n","processBismark timing...")
  end_time <- Sys.time()
  print(end_time - start_time)
  
  print(glue::glue("Before filtering there were {nrow(bs)} CpGs, \\
                  after filtering for {Cov}x coverage in at least {per.ctrl*100} ctrl samples and {per.exp*100} exp samples, \\
                  there are now {nrow(bs.filtered)} CpGs assayed"))
  
  return(bs.filtered)
}
