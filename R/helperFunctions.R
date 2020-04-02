#' getSmooth
#' @description Provides individual smoothed methylation values for a \code{GRanges} object using \code{bsseq}
#' @param bsseq A smoothed \code{bsseq} object
#' @param regions A \code{GRanges} object of regions to obtain smoothed methylation values for
#' @return A data frame of individual smoothed methylation values
#' @import bsseq
#' @importFrom glue glue
#' @export getSmooth
getSmooth <- function(bsseq = bs.filtered.bsseq,
                      regions = sigRegions){
  print(glue::glue("Obtaining smoothed methylation values..."))
  data.frame(
    bsseq::getMeth(BSseq = bsseq,
                   regions = regions,
                   type = "smooth",
                   what = "perRegion"),
    check.names = FALSE) %>% 
    cbind(regions, .) %>% 
    return()
}

#' smooth2txt
#' @description Save smoothed methylation values as a text file
#' @param df Data frame
#' @param txt Character string of save file name
#' @return Saves a text file
#' @importFrom glue glue
#' @export smooth2txt
smooth2txt <- function(df = df,
                       txt = txt){
  print(glue::glue("Saving individual smoothed methylation values to {txt}"))
  write.table(df,
              txt,
              sep = "\t",
              quote = FALSE,
              row.names = FALSE,
              col.names = TRUE)
}

#' gr2csv
#' @description Save a genomic ranges object as a csv file
#' @param gr \code{GRanges} or \code{bsseq} object
#' @param csv Character string of save file name
#' @return Saves a CSV file
#' @import BiocGenerics
#' @importFrom glue glue
#' @export gr2csv
gr2csv <- function(gr = gr,
                   csv = csv){
  print(glue::glue("Saving {csv}"))
  write.csv(BiocGenerics::as.data.frame(gr),
            file = csv,
            row.names = FALSE)
}

#' gr2bed
#' @description Save a genomic ranges object as a basic bed file
#' @param gr Genomic ranges or bsseq object
#' @param bed Name of the bed file in quotations
#' @return Bed file
#' @import BiocGenerics
#' @importFrom glue glue
#' @export gr2bed
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

#' df2bed
#' @description Save a dataframe as a basic bed file
#' @param df Data frame
#' @param bed Name of the bed file in quotations
#' @return Bed file
#' @importFrom glue glue
#' @export df2bed
df2bed <-function(df = df,
                  bed = bed){
  print(glue::glue("Saving {bed}"))
  write.table(df,
              bed,
              quote = FALSE,
              row.names = FALSE,
              col.names = FALSE,
              sep = "\t")
}

#' gg_color_hue
#' @description Generate ggplot2 style colors
#' @param n Number of samples
#' @return Character string of colors
#' @references \url{https://stackoverflow.com/questions/8197559/emulate-ggplot2-default-color-palette}
#' @importFrom glue glue
#' @export gg_color_hue
gg_color_hue <- function(n = n){
  print(glue::glue("Preparing colors for {n} samples"))
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

#' getBackground
#' @description Get background regions from filtered bsseq object based on minCpGs and maxGap
#' @param bs bsseq object that has been filtered for coverage and sorted
#' @param minNumRegion Minimum CpGs required for a region. Must be at least 3 CpGs, default is 5 CpGs
#' @param maxGap Maximum distance between CpGs to be included in the same region. Default is 1000 bp
#' @return Data.frame of background regions with location, number of CpGs, and width
#' @import bsseq
#' @export getBackground
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

#' manQQ
#' @description Create manhattan and Quantile-Quantile (Q-Q) plots of \code{ChIPseeker csAnno} peak object with genic annotations using \code{CMplot}
#' @param peakAnno A \code{ChIPseeker csAnno} peak object of background regions from \code{dmrseq::dmrseq()}
#' @param ... Additional arguments passed onto \code{\link[CMplot]{CMplot}}
#' @return Saves a pdf of manhattan and qq plots
#' @import CMplot
#' @import ChIPseeker
#' @import GenomicRanges
#' @import RColorBrewer
#' @importFrom glue glue
#' @export manQQ
manQQ <- function(backgroundAnno = backgroundAnno,
                  ...){
  cat("\n[DMRichR] Manhattan and QQ plots \t\t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
  glue::glue("Generating Manhattan and QQ plots...")
  setwd("DMRs")
  CMplot::CMplot(backgroundAnno %>%
                   ChIPseeker::as.GRanges() %>%
                   sort() %>%
                   as.data.frame() %>%
                   dplyr::select(SYMBOL, seqnames, start, pval) %>%
                   dplyr::mutate(seqnames = substring(.$seqnames, 4)),
                 col = DMRichR::gg_color_hue(2),
                 plot.type = c("m","q"),
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
#' @description Save DMRs and background regions from \code{dmrseq::dmrseq()} in formats for external analyses using GAT and HOMER
#' @param siRegions A \code{GRanges} object of signficant DMRs returned by \code{dmrseq::dmrseq()}
#' @param regions A \code{GRanges} object of background regions returned by \code{dmrseq::dmrseq()}
#' @return Saves external GAT and HOMER folders with bed files into the working directory
#' @import tidyverse
#' @importFrom glue glue
#' @export saveExternal
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

#' lift450k
#' @description LiftOver 450k CpG IDs to hg38 coordinates
#' @param probes A dataframe or vector of 450k CpG IDss
#' @return A \code{GRanges} object of hg38 coordinates
#' @import tidyverse
#' @import FDb.InfiniumMethylation.hg19
#' @import rtracklayer
#' @import R.utils
#' @import GenomicRanges
#' @importFrom glue glue
#' @export lift450k
lift450k <- function(probes = probes){
  
  hm450 <- get450k()
  
  if(!file.exists("hg19ToHg38.over.chain")){
    message("Downloading hg19 to hg38 LiftOver chain..")
    url <- "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz"
    download.file(url, basename(url))
    gunzip(basename(url))
  }
  
  message(glue::glue("Performing liftOver to hg38 for {nrow(probes)} probes..."))
  
  if(is.data.frame(probes)){
    probes <- probes %>%
      dplyr::pull()
  }
  
  hg19 <- hm450[probes][,0]
  chain <- import.chain("hg19ToHg38.over.chain")
  hg38 <- unlist(liftOver(hg19, chain))
  
  message(glue::glue("{length(hg38)} out of {length(hg19)} probes were liftedOver..."))
  
  return(hg38)
}
