#' getBackground
#' @title bsseq style background regions
#' @description Get background regions from filtered bsseq object based on minCpGs and maxGap
#' @param bs bsseq object that has been filtered for coverage and sorted
#' @param minNumRegion Minimum CpGs required for a region. Must be at least 3 CpGs, default is 5 CpGs
#' @param maxGap Maximum distance between CpGs to be included in the same region. Default is 1000 bp
#' @return Data.frame of background regions with location, number of CpGs, and width
#' @import bsseq
#' @export getBackground
#' 
getBackground <- function(bs = bs.filtered,
                          minNumRegion = 5,
                          maxGap = 1000){
        background <- bsseq:::regionFinder3(x = as.integer(rep(1, length(bs))), chr = as.character(seqnames(bs)), 
                                            positions = start(bs), maxGap = maxGap, verbose = FALSE)[["up"]]
        background <- subset(background, n >= minNumRegion, select = c("chr", "start", "end", "n"))
        background$chr <- as.character(background$chr)
        background$width <- background$end - background$start
        return(background)
}

#' Manhattan
#' @title Manhattan plot
#' @description Create a Manhattan plot of \code{ChIPseeker csAnno} peak object with genic annotations using \code{CMplot}
#' @param backgroundAnno A \code{ChIPseeker csAnno} peak object of background regions from \code{dmrseq::dmrseq()}
#' @param ... Additional arguments passed onto \code{\link[CMplot]{CMplot}}
#' @return Saves a pdf of manhattan and qq plots
#' @import CMplot
#' @import ChIPseeker
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom RColorBrewer brewer.pal
#' @importFrom glue glue
#' @importFrom magrittr %>%
#' @importFrom dplyr as_tibble select mutate 
#' @export Manhattan
#' 
Manhattan<- function(backgroundAnno = backgroundAnno,
                     ...){
  cat("\n[DMRichR] Manhattan plot \t\t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
  setwd("DMRs")
  CMplot::CMplot(backgroundAnno %>%
                   GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE) %>%
                   sort() %>%
                   dplyr::as_tibble() %>%
                   dplyr::select(geneSymbol, seqnames, start, p.value) %>%
                   dplyr::mutate(seqnames = substring(.$seqnames, 4)),
                 col = DMRichR::gg_color_hue(2),
                 plot.type = "m",
                 LOG10 = TRUE,
                 ylim = NULL,
                 threshold = 0.05, #c(1e-6,1e-4),
                 threshold.lty = c(1,2),
                 threshold.lwd = c(1,1),
                 threshold.col = c("black","grey"),
                 cex = 0.5,
                 cex.axis = 0.7,
                 amplify = FALSE,
                 chr.den.col = RColorBrewer::brewer.pal(9, "YlOrRd"),
                 bin.size = 1e6,
                 signal.col = c("red","green"),
                 signal.cex = c(1,1),
                 signal.pch = c(19,19),
                 file = "pdf",
                 memo = "",
                 ...)
  setwd('..')
}

#' saveExternal
#' @title Save regions for external enrichment testing
#' @description Save DMRs and background regions from \code{dmrseq::dmrseq()} in formats for external analyses using GAT and HOMER
#' @param sigRegions A \code{GRanges} object of signficant DMRs returned by \code{dmrseq::dmrseq()}
#' @param regions A \code{GRanges} object of background regions returned by \code{dmrseq::dmrseq()}
#' @return Saves external GAT and HOMER folders with bed files into the working directory
#' @importFrom glue glue
#' @importFrom magrittr %>%
#' @importFrom dplyr as_tibble mutate case_when select filter
#' @export saveExternal
#' 
saveExternal <- function(sigRegions = sigRegions,
                         regions = regions){
  cat("\n[DMRichR] Preparing files for annotations \t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
  
  if(dir.exists("Extra") == F){dir.create("Extra")}
  
  glue::glue("Preparing regions for external GAT analysis...")
  dir.create("Extra/GAT")
  
  sigRegions <- sigRegions %>% 
    dplyr::as_tibble() %>%
    dplyr::mutate(direction = dplyr::case_when(stat > 0 ~ "Hypermethylated",
                                               stat < 0 ~ "Hypomethylated"
                                               )
                  )
  
  regions <- regions %>%
    dplyr::as_tibble() 
  
  sigRegions %>%
    dplyr::select(seqnames, start, end, direction) %>% 
    DMRichR::df2bed("Extra/GAT/DMRs.bed")
  
  regions %>%
    dplyr::select(seqnames, start, end) %>% 
    DMRichR::df2bed("Extra/GAT/background.bed")
  
  glue::glue("Preparing DMRs for external HOMER analysis...")
  dir.create("Extra/HOMER")
  
  sigRegions %>%
    dplyr::filter(direction == "Hypermethylated") %>%
    dplyr::select(seqnames, start, end) %>%
    DMRichR::df2bed("Extra/HOMER/DMRs_hyper.bed")
  
  sigRegions %>%
    dplyr::filter(direction == "Hypomethylated") %>%
    dplyr::select(seqnames, start, end) %>%
    DMRichR::df2bed("Extra/HOMER/DMRs_hypo.bed")
  
  regions %>%
    dplyr::select(seqnames, start, end) %>% 
    DMRichR::df2bed("Extra/HOMER/background.bed")
  
}

#' arrayLift
#' @title LiftOver Infinium array probe IDs to hg38 coordinates
#' @description LiftOver EPIC, 450k, or 27K Infinium array CpG IDs to hg38 coordinates
#' @param probes A dataframe or vector of EPIC, 450K, or 27K CpG IDs
#' @param array A character with array platform ("EPIC", "450K" or "27K")
#' @return A \code{GRanges} object of hg38 coordinates
#' @import GenomicRanges
#' @importFrom dplyr as_tibble select pull
#' @importFrom tibble rownames_to_column
#' @importFrom rtracklayer liftOver
#' @importFrom AnnotationHub AnnotationHub
#' @importFrom glue glue
#' @importFrom magrittr %>%
#' @importFrom minfi getAnnotation
#' @importFrom plyranges mutate filter
#' @references \url{https://support.bioconductor.org/p/78652/}
#' @examples
#' \dontrun{ 
#' readxl::read_excel("file.xlsx") %>%
#'   dplyr::select(CpGID) %>% 
#'   arrayLift("EPIC") %>%
#'   dplyr::as_tibble() %>%
#'   dplyr::select(seqnames, start, end) %>% 
#'   DMRichR::df2bed("file.bed")
#'}
#' @export arrayLift
#' 
arrayLift <- function(probes = probes,
                      array = "EPIC"){
  glue::glue("Obtaining probes from {array}")
  if(array == "EPIC"){
    
    if(!require(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)){
      BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")}
    library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
    
    message("Fetching coordinates for hg19...")
    array <- minfi::getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19) %>%
      as.data.frame() %>% 
      tibble::rownames_to_column() %>% 
      dplyr::select(rowname, chr, pos, strand) %>% 
      GenomicRanges::makeGRangesFromDataFrame(.,seqnames.field = "chr",
                                              start.field = "pos",
                                              end.field = "pos",
                                              strand.field = "strand",
                                              keep.extra.columns = TRUE) %>% 
      DMRichR::extend(downstream = 1)
    
  }else if(array == "450K"){
    
    if(!require(FDb.InfiniumMethylation.hg19)){
      BiocManager::install("FDb.InfiniumMethylation.hg19")}
    library(FDb.InfiniumMethylation.hg19)
    
    array <- FDb.InfiniumMethylation.hg19::get450k() %>%
      plyranges::mutate(rowname = names(.))
    
  }else if(array == "27K"){
    
    if(!require(FDb.InfiniumMethylation.hg19)){
      BiocManager::install("FDb.InfiniumMethylation.hg19")}
    library(FDb.InfiniumMethylation.hg19)
    
    array <- FDb.InfiniumMethylation.hg19::get27k() %>%
      plyranges::mutate(rowname = names(.))
    
  }else{
    stop(glue::glue("{array} is not suppourted, please choose EPIC, 450K, or 27K"))
  }
  
  message("Performing liftOver to hg38...")
  
  if(is.data.frame(probes)){
    probes <- probes %>%
      dplyr::pull()
  }
  
  hg38 <- array %>%
    plyranges::filter(rowname %in% probes) %>% 
    rtracklayer::liftOver(AnnotationHub::AnnotationHub()[["AH14150"]]) %>%
    unlist()
  
  message(glue::glue("{length(hg38)} out of {length(probes)} probes were liftedOver..."))
  
  return(hg38)
}

#' extend
#' @title Extend genomic ranges
#' @description Extend the 5' and 3' ends of a \code{GRanges}. See references for source. 
#' @param x A \code{GRanges} object to extend
#' @param upstream Numeric of basepairs to extend the 5' end by
#' @param downstream Numeric of basepairs to extend the 3' end by
#' @return A \code{GRanges} object with extended ranges
#' @import GenomicRanges
#' @references \url{https://support.bioconductor.org/p/78652/}
#' @export extend
#' 
extend <- function(x,
                   upstream = 0,
                   downstream = 0)
{
  if (any(strand(x) == "*"))
    warning("'*' ranges were treated as '+'")
  on_plus <- strand(x) == "+" | strand(x) == "*"
  new_start <- start(x) - ifelse(on_plus, upstream, downstream)
  new_end <- end(x) + ifelse(on_plus, downstream, upstream)
  ranges(x) <- IRanges(new_start, new_end)
  trim(x)
}

#' read_excel_all
#' @title Read entire excel document
#' @description Read all sheets in an excel document 
#' @param filename A character vector specifying the name of the excel document
#' @return A list of \code{tibbles} with the excel document data 
#' @importFrom readxl excel_sheets read_excel
#' @importFrom purrr set_names map 
#' @references \url{https://stackoverflow.com/a/12945838}
#' @export read_excel_all
#' 
read_excel_all <- function(filename) {
  readxl::excel_sheets(filename) %>% 
    purrr::set_names() %>%
    purrr::map(function(x){readxl::read_excel(filename, sheet = x)}) %>%
    return()
}

#' smooth2txt
#' @title Save regions and methylation values
#' @description Provides individual smoothed methylation values for a \code{GRanges} object using \code{bsseq}
#' @param bsseq A smoothed \code{bsseq} object
#' @param regions A \code{GRanges} object of regions to obtain smoothed methylation values for
#' @param txt Character string of save file name
#' @return Saves a text file
#' @importFrom bsseq getMeth
#' @importClassesFrom bsseq BSseq 
#' @importFrom magrittr %>%
#' @importFrom glue glue
#' @importFrom utils write.table
#' @export smooth2txt
#' 
smooth2txt <- function(bsseq = bs.filtered.bsseq,
                       regions = sigRegions,
                       txt = txt){
  print(glue::glue("Saving individual smoothed methylation values to {txt}"))
  data.frame(
    bsseq::getMeth(BSseq = bsseq,
                   regions = regions,
                   type = "smooth",
                   what = "perRegion"),
    check.names = FALSE) %>% 
    cbind(regions, .) %>% 
    write.table(.,
                txt,
                sep = "\t",
                quote = FALSE,
                row.names = FALSE,
                col.names = TRUE)
}

#' gr2bed
#' @title Save a genomic ranges object as a bed file
#' @description Save a genomic ranges object as a basic bed file
#' @param gr Genomic ranges or bsseq object
#' @param bed Name of the bed file in quotations
#' @return Bed file
#' @importFrom BiocGenerics as.data.frame
#' @importFrom glue glue
#' @export gr2bed
#' 
gr2bed <- function(gr = gr,
                   bed = bed){
  print(glue::glue("Saving {bed}"))
  write.table(BiocGenerics::as.data.frame(gr)[1:3],
              bed,
              sep = "\t",
              row.names = FALSE,
              col.names = FALSE,
              quote = FALSE)
}

#' gg_color_hue
#' @title ggplot2 colors
#' @description Generate ggplot2 style colors
#' @param n Number of samples
#' @return Character string of colors
#' @references \url{https://stackoverflow.com/questions/8197559/emulate-ggplot2-default-color-palette}
#' @importFrom glue glue
#' @export gg_color_hue
#' 
gg_color_hue <- function(n = n){
  print(glue::glue("Preparing colors for {n} samples"))
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

