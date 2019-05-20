#' GREAT
#' @description Perform gene ontology of DMRs relative to background regions from \code{dmrseq::dmrseq()} using \code{rGREAT}
#' @param siRegions A \code{GRanges} object of signficant DMRs returned by \code{dmrseq:dmrseq()}
#' @param regions A \code{GRanges} object of background regions returned by \code{dmrseq:dmrseq()}
#' @param genome A character vector specifying the genome of interest ("hg38" or "mm10")
#' @return A list of gene ontology enrichments.
#' @import rGREAT
#' @import liftOver
#' @import rtracklayer
#' @export GREAT
GREAT <- function(sigRegions = sigRegions,
                  regions = regions, 
                  genome = genome){
  if(genome == "hg38"){
    cat("\n[DMRichR] Obtaining liftOver information for GREAT \t", format(Sys.time(), "%d-%m-%Y %X"))
    path <- system.file(package="liftOver", "extdata", "hg38ToHg19.over.chain")
    ch <- import.chain(path)
    ch
    
    glue::glue("liftOver DMRs...")
    seqlevelsStyle(sigRegions) <- "UCSC"
    sigRegions_liftOver <- liftOver(sigRegions, ch)
    class(sigRegions_liftOver)
    sigRegions_liftOver <- unlist(sigRegions_liftOver)
    #stopifnot(length(sigRegions) > length(sigRegions_liftOver))
    
    glue::glue("liftOver background regions...")
    seqlevelsStyle(regions) <- "UCSC"
    regions_liftOver <- liftOver(regions, ch)
    class(regions_liftOver)
    regions_liftOver <- unlist(regions_liftOver)
    #stopifnot(length(regions) > length(regions_liftOver))
  }
  
  if(genome == "hg38" | genome == "mm10"){
    cat("\n[DMRichR] Submitting to GREAT\ \t\t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
    if(genome == "hg38"){
      gr <- sigRegions_liftOver
      bg <- regions_liftOver
      species <- "hg19"
    }else if(genome == "mm10"){
      gr <- sigRegions
      bg <- regions
      species <- "mm10"
    }else{
      stop(paste(genome, "is not suppourted for GREAT, please choose either hg38 or mm10 [Case Sensitive]"))
    }
    
    job <- submitGreatJob(gr,
                          bg = bg,
                          species = species,
                          request_interval = 1)
    job
    #availableCategories(job)
    #availableOntologies(job)
  }
  return(job)
}
