#' prepareHOMER
#' @title Save regions for HOMER
#' @description Save DMRs and background regions from \code{dmrseq::dmrseq()}for HOMER
#' @param sigRegions A \code{GRanges} object of significant DMRs returned by \code{dmrseq::dmrseq()}
#' @param regions A \code{GRanges} object of background regions returned by \code{dmrseq::dmrseq()}
#' @return Creates a folder for HOMER with bed files
#' @importFrom magrittr %>%
#' @importFrom plyranges filter
#' @export prepareHOMER
#' 
prepareHOMER <- function(sigRegions = sigRegions,
                         regions = regions){
  
  dir.create("HOMER")
  
  sigRegions %>%
    DMRichR::gr2bed("HOMER/DMRs.bed")
  
  sigRegions %>%
    plyranges::filter(stat > 0) %>% 
    DMRichR::gr2bed("HOMER/DMRs_hyper.bed")
  
  sigRegions %>%
    plyranges::filter(stat < 0) %>% 
    DMRichR::gr2bed("HOMER/DMRs_hypo.bed")
  
  regions %>%
    DMRichR::gr2bed("HOMER/background.bed")
  
}

#' HOMER
#' @title Transcription factor motif analysis
#' @description Run a HOMER known transcription motif analysis for all DMRs,
#' hypermethylated DMRs, and hypomethylated DMRs. The function requires HOMER's 
#' \code{findMotifsGenome.pl} script to be in the path (i.e. \code{module load homer})
#' and for the genome of interest (see \code{perl /path-to-homer/configureHomer.pl -list}) 
#' to be installed through \code{perl /path-to-homer/configureHomer.pl -install human}.
#' @param genome Character specifying the genome
#' @param cores Integer specifying the number of cores to use
#' @return A folder with HOMER results
#' @importFrom glue glue
#' @references \url{http://homer.ucsd.edu/homer/introduction/configure.html}
#' @references \url{https://www.biostars.org/p/443759/}
#' @export HOMER
#' 
HOMER <- function(genome = genome,
                  cores = cores){
  cat("\n[DMRichR] HOMER known transcription factor motif analysis \t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
  tryCatch({
    if(Sys.which("findMotifsGenome.pl") == ""){
      print(glue::glue("HOMER was not detected in PATH, skipping motif analysis. Did you load the module?"))
    }else{
      print(glue::glue("HOMER was detected in PATH, now performing motif enrichment for {genome} using {cores} cores"))
      system(paste(shQuote(system.file("exec/HOMER.sh", package = "DMRichR")),genome,cores))
    }
  },
  error = function(error_condition) {
    print(glue::glue("There was an error with running HOMER for {genome}. 
                    Have you confirmed that {genome} is avaiable in HOMER and installed using:
                    perl {path}/configureHomer.pl -list
                    perl {path}/configureHomer.pl -install {genome}
                    Note: see https://www.biostars.org/p/443759/ for HOMER path if using a conda install.",
                     path = dirname(Sys.which("findMotifsGenome.pl"))))
  })
}
