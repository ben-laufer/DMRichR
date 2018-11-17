#' packageManage
#' @description Install package management
#' @export packageManage
packageManage <- function(){
  CRAN <- c("BiocManager", "remotes")
  new.packages <- CRAN[!(CRAN %in% installed.packages()[,"Package"])]
  if(length(new.packages)>0){
    install.packages(new.packages, repos ="https://cloud.r-project.org", quiet = TRUE)
  }
  stopifnot(suppressMessages(sapply(CRAN, require, character.only = TRUE)))
}

#' packageLoad
#' @description Install and load desired packages
#' @param packages Character string of desired packages
#' @export packageLoad
packageLoad <- function(packages = packages){
  new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)>0){
    library("BiocManager")
    new.packages <- gsub("ggbiplot", "vqv/ggbiplot", packages)
    BiocManager::install(new.packages, ask = FALSE, quiet = TRUE)
  }
  stopifnot(suppressMessages(sapply(packages, require, character.only = TRUE)))
}

#' getSmooth
#' @description Provides individual smoothed methylation values for genomic ranges objects using bsseq
#' @param bsseq Smoothed bsseq object
#' @param regions Genomic ranges object
#' @param out Name of the text file in quotations
#' @return Genomic ranges object of individual smoothed methylation values and text file
#' @export getSmooth
getSmooth <- function(bsseq = bsseq,
                      regions = regions,
                      out = out){
  message("Smoothing...")
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
  message("Saving bed file...")
  write.table(df, txt, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}


#' gr2csv
#' @description Save a genomic ranges object as a csv file
#' @param gr Genomic ranges or bsseq object
#' @param csv Name of the csv file in quotations
#' @return CSV file
#' @export gr2csv
gr2csv <- function(gr = gr,
                   csv = csv){
  message("Saving CSV...")
  write.csv(as.data.frame(gr), file = csv, row.names = FALSE)
}

#' gr2bed
#' @description Save a genomic ranges object as a basic bed file
#' @param gr Genomic ranges or bsseq object
#' @param bed Name of the bed file in quotations
#' @return Bed file
#' @export gr2bed
gr2bed <- function(gr = gr,
                   bed = bed){
  message("Saving bed file...")
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
  message("Saving bed file...")
  write.table(df, bed, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
}

#' gg_color_hue
#' @description Generate ggplot2 style colors
#' @param n Number of samples
#' @return Character string of colors
#' @references \url{https://stackoverflow.com/questions/8197559/emulate-ggplot2-default-color-palette}
#' @export gg_color_hue
gg_color_hue <- function(n = n){
  message("Preparing colors...")
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
