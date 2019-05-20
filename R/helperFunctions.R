#' packageLoad
#' @description Install and load desired packages
#' @param packages Character string of desired packages
#' @export packageLoad
packageLoad <- function(packages = packages){
  cat("\n","Checking for BiocManager and helpers...")
  CRAN <- c("BiocManager", "remotes", "magrittr")
  new.CRAN.packages <- CRAN[!(CRAN %in% installed.packages()[,"Package"])]
  if(length(new.CRAN.packages)>0){
    install.packages(new.CRAN.packages, repos ="https://cloud.r-project.org", quiet = TRUE)
  }
  cat("Done")
  cat("\n", "Loading package management...")
  stopifnot(suppressMessages(sapply(CRAN, require, character.only = TRUE)))
  cat("Done")
  
  new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)>0){
    cat("Installing missing packages...")
    new.packages <- packages %>%
      gsub("ggbiplot", "vqv/ggbiplot", .) %>% 
      gsub("DMRichR", "ben-laufer/DMRichR", .) %>% 
      gsub("gt", "rstudio/gt", .)
    BiocManager::install(new.packages, ask = FALSE, quiet = TRUE)
    cat("Done")
  }
  cat("\n", "Loading packages...")
  stopifnot(suppressMessages(sapply(packages, require, character.only = TRUE)))
  suppressWarnings(BiocManager::valid(fix = TRUE, update = TRUE, ask = FALSE))
  cat("Done", "\n")
}

#' getSmooth
#' @description Provides individual smoothed methylation values for genomic ranges objects using bsseq
#' @param bsseq Smoothed bsseq object
#' @param regions Genomic ranges object
#' @param out Name of the text file in quotations
#' @return Genomic ranges object of individual smoothed methylation values and text file
#' @import bsseq
#' @export getSmooth
getSmooth <- function(bsseq = bsseq,
                      regions = regions,
                      out = out){
  print(glue::glue("Obtaining smoothed methylation values..."))
  smoothed <- data.frame(getMeth(BSseq = bsseq, regions = regions, type = "smooth", what = "perRegion"), check.names=FALSE)
  smoothed_table <- cbind(regions, smoothed)
  write.table(smoothed_table, out, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  return(smoothed_table)
}

#' smooth2txt
#' @description Save smoothed methylation values as a text file
#' @param df Data frame
#' @param txt Name of the text file in quotations
#' @return Text file
#' @export smooth2txt
smooth2txt <- function(df = df,
                       txt = txt){
  print(glue::glue("Saving bed file..."))
  write.table(df, txt, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}


#' gr2csv
#' @description Save a genomic ranges object as a csv file
#' @param gr Genomic ranges or bsseq object
#' @param csv Name of the csv file in quotations
#' @return CSV file
#' @import BiocGenerics
#' @export gr2csv
gr2csv <- function(gr = gr,
                   csv = csv){
  print(glue::glue("Saving CSV..."))
  write.csv(as.data.frame(gr), file = csv, row.names = FALSE)
}

#' gr2bed
#' @description Save a genomic ranges object as a basic bed file
#' @param gr Genomic ranges or bsseq object
#' @param bed Name of the bed file in quotations
#' @return Bed file
#' @import BiocGenerics
#' @export gr2bed
gr2bed <- function(gr = gr,
                   bed = bed){
  print(glue::glue("Saving bed file..."))
  write.table(as.data.frame(gr)[1:3], bed, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
}

#' df2bed
#' @description Save a dataframe as a basic bed file
#' @param df Data frame
#' @param bed Name of the bed file in quotations
#' @return Bed file
#' @export df2bed
df2bed <-function(df = df,
                  bed = bed){
  print(glue::glue("Saving bed file..."))
  write.table(df, bed, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
}

#' gg_color_hue
#' @description Generate ggplot2 style colors
#' @param n Number of samples
#' @return Character string of colors
#' @references \url{https://stackoverflow.com/questions/8197559/emulate-ggplot2-default-color-palette}
#' @export gg_color_hue
gg_color_hue <- function(n = n){
  print(glue::glue("Preparing colors..."))
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
        background$start <- background$start - 1
        background$end <- background$end + 1
        background$width <- background$end - background$start
        return(background)
}

#' labelDirection
#' @description Annotate a \code{GRanges} object containing either DMRs or background regions with directionality of change
#'  (i.e. hypermethlated or hypomethylated).
#' @param regions A \code{GRanges} object returned by \code{dmrseq:dmrseq()}.
#' @return An \code{GRanges} object annotated with directionality of change.
#' @import GenomicRanges
#' @export labelDirection
labelDirection <- function(regions = sigRegions){
  cat("\n", "Annotating regions with directionality...")
  external <- regions
  for (i in 1:length(external)){
    if(external$stat[i] > 0){
      external$direction[i] <- "Hypermethylated"
    }else if(external$stat[i] < 0){
      external$direction[i] <- "Hypomethylated"
    }else{
      stop("Annotation problem")
    }}
  cat("Done", "\n")
  return(external)
}

#' tidyDMRs
#' @description Tidy DMRs or background regions that have been annotated using ChIPseeker
#' @param regions Peak file or GRanges object from ChIPseeker
#' @import ChIPseeker
#' @import tidyverse
#' @export tidyDMRs
tidyDMRs <- function(regions = peakAnno){
  regions %>% 
    as.tibble() %>%
    dplyr::select("seqnames", "start", "end", "width", "L",
                  "beta", "stat", "pval", "qval", "percentDifference",
                  "annotation", "distanceToTSS", "ENSEMBL", "SYMBOL", "GENENAME") %>%
    dplyr::rename(CpGs = L,
                  betaCoefficient = beta, statistic = stat,
                  "p-value" = pval, "q-value" = qval, difference = percentDifference,
                  geneSymbol = SYMBOL, gene = GENENAME) %>% 
    return()
}

#' DMReport
#' @description Create an html report of a \code{ChIPseeker csAnno} peak object with genic annotations.
#' @param peakAnno A \code{ChIPseeker csAnno} peak object of DMRs from \code{dmrseq::dmrseq()}.
#' @return Saves an html report of DMRs with genic annotations.
#' @import gt
#' @export DMReport
DMReport <- function(peakAnno = peakAnno){
  cat("\n","Preparing HTML report...")
  peakAnno %>%
    DMRichR::tidyDMRs() %>% 
    dplyr::select(-ENSEMBL, -betaCoefficient, -statistic) %>%
    dplyr::mutate(difference = difference/100) %>% 
    gt() %>%
    tab_header(
      title = glue::glue("{length(sigRegions)} Significant DMRs"),
      subtitle = glue::glue("{length(sigRegions)} Significant DMRs \\
                         ({round(sum(sigRegions$stat > 0) / length(sigRegions), digits = 2)*100}% hypermethylated, \\
                         {round(sum(sigRegions$stat < 0) / length(sigRegions), digits = 2)*100}% hypomethylated) \\
                         in {length(regions)} background regions \\
                         from {nrow(bs.filtered)} CpGs assayed at {coverage}x coverage")
    ) %>% 
    fmt_number(
      columns = vars("width", "CpGs"),
      decimals = 0
    ) %>% 
    fmt_scientific(
      columns = vars("p-value", "q-value"),
      decimals = 2
    ) %>%
    fmt_percent(
      columns = vars("difference"),
      drop_trailing_zeros = TRUE
    ) %>% 
    as_raw_html(inline_css = TRUE) %>%
    write("DMRs.html")
  cat("Done", "\n")
}

#' manQQ
#' @description Create manhattan and Quantile-Quantile (Q-Q) plots of \code{ChIPseeker csAnno} peak object with genic annotations using \code{CMplot}
#' @param peakAnno A \code{ChIPseeker csAnno} peak object of background regions from \code{dmrseq::dmrseq()}
#' @param ... Additional arguments passed onto \code{CMplot::CMplot}
#' @return Saves a pdf of manhattan and qq plots
#' @import CMplot
#' @import ChIPseeker
#' @import GenomicRanges
#' @export manQQ
manQQ <- function(backgroundAnno = backgroundAnno,
                  ...){
  cat("\n[DMRichR] Manhattan and QQ plots \t\t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
  glue::glue("Tidying for Manhattan and QQ plots")
  Manhattan <- as.data.frame(sort(as.GRanges(backgroundAnno), ignore.strand=TRUE))[c("SYMBOL","seqnames", "start", "pval")]
  Manhattan$seqnames <- substring(Manhattan$seqnames, 4)
  
  glue::glue("Generating Manhattan and QQ plots...")
  CMplot(Manhattan,
         col = gg_color_hue(2),
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
         chr.den.col = brewer.pal(9, "YlOrRd"),
         bin.size = 1e6,
         bin.max = 100,
         signal.col = c("red","green"),
         signal.cex = c(1,1),
         signal.pch = c(19,19),
         file = "pdf",
         memo = "",
         ...)
}
