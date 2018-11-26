#' smoothHeatmap
#' @description Plot a heatmap of normalized individual smoothed methylation value z scores for selected regions (i.e. significant DMRs)
#' @param bsseq Smoothed bsseq object
#' @param regions Genomic ranges object
#' @param names Ordered sample names
#' @param groups Ordered test covariate information for each sample
#' @param out Name of the text file to save in quotations
#' @param ... Additional arguments passed onto heatmap.3()
#' @return Saves a pdf image of the heatmap
#' @import bsseq
#' @import tidyverse
#' @import gplots
#' @references \url{https://sebastianraschka.com/Articles/heatmaps_in_r.html}
#' @references \url{https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R}
#' @references \url{https://www.biostars.org/p/18211/}
#' @export smoothHeatmap
smoothHeatmap <- function(regions = sigRegions,
                          bsseq = bs.filtered.bsseq,
                          groups = bs.filtered.bsseq %>% pData() %>% as.tibble() %>% pull(!!testCovariate),
                          out = "sig_individual_smoothed_DMR_methylation.txt",
                          ...){
  cat("\n[DMRichR] DMR heatmap \t\t\t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
  pDataFactors <- pData(bsseq) %>% as.data.frame() %>% dplyr::select_if(is.factor)
  ColSideColors <- matrix(nrow = nrow(pDataFactors), ncol = ncol(pDataFactors))
  
  for(i in 1:length(pDataFactors)){
    matrix <- c(gg_color_hue(length(levels(pDataFactors[,i]))))[pDataFactors[,i]]
    matrix  <- mapvalues(pDataFactors[,i],
                         from = levels(pDataFactors[,i]),
                         to = unique(matrix)) %>% 
      as.matrix()
    ColSideColors[,i] <- matrix 
  }
  colnames(ColSideColors) <- names(pDataFactors)
  
  
  message("Obtaining smoothed methylation values...")
  smoothed <- data.frame(getMeth(BSseq = bsseq, regions = regions, type = "smooth", what = "perRegion"))
  smoothed_table <- cbind(regions, smoothed)
  write.table(smoothed_table, out, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  
  message("Tidying for heatmap of HCA...")
  # Load smoothed values
  matrix <- as.matrix(smoothed)
  # Convert to Percent
  matrix <- matrix[,]*100
  # Subtract the mean methylation for each row/DMR
  data <- sweep(matrix, 1, rowMeans(matrix))
  # Tidy
  data <- as.matrix(data)
  colnames(data) <- groups
  
  message("Plotting heatmap of HCA...")
  source("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")
  heatmap.3(data,
            Rowv= as.dendrogram(hclust(dist(data))),
            scale = c("row"),
            Colv = TRUE,
            col = rev(brewer.pal(11, name = "RdBu")),
            margins = c(10,10),
            trace = "none",
            main = paste(length(regions),"Differentially Methylated Regions", sep = " "),
            labRow = NA,
            #keysize = 0.85,
            #key.par = list(cex=0.5),
            KeyValueName = "Z-score(% mCG/CG - mean)",
            ColSideColors = ColSideColors,
            ...
  )
  
  # par(lend = 1)
  # legend("topright",
  #        legend = levels(pDataFactor),
  #        col = levels(ColSideColors), 
  #        lty= 1, 
  #        lwd = 10)
  
}
