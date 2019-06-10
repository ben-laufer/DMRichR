#' smoothHeatmap
#' @description Plot a heatmap of normalized individual smoothed methylation value z scores for selected regions (i.e. significant DMRs)
#' The plotted values are not typical Z-scores, but rather a Z-score of %mCG/CG - mean,
#'  which is data visualization trick that allows you to focus on the variance
#' @param bsseq Smoothed \code{bsseq} object
#' @param regions \code{GRanges} object of regions to plot a heatmap for
#' @param groups Ordered test covariate information for each sample
#' @param ... Additional arguments passed onto \code{heatmap.3()}
#' @return An image of the heatmap
#' @import bsseq
#' @import gplots
#' @importFrom plyr mapvalues
#' @references \url{https://sebastianraschka.com/Articles/heatmaps_in_r.html}
#' @references \url{https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R}
#' @references \url{https://www.biostars.org/p/18211/}
#' @export smoothHeatmap
smoothHeatmap <- function(regions = sigRegions,
                          bsseq = bs.filtered.bsseq,
                          groups = bs.filtered.bsseq %>% pData() %>% dplyr::as_tibble() %>% dplyr::pull(!!testCovariate),
                          ...){
  cat("\n[DMRichR] DMR heatmap \t\t\t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
  
  print(glue::glue("Determining factor colors..."))
  pDataFactors <- pData(bsseq) %>%
    as.data.frame() %>%
    dplyr::select_if(is.factor)
  ColSideColors <- matrix(nrow = nrow(pDataFactors), ncol = ncol(pDataFactors))
  legendNames <- as.character()
  legendColors <- as.character()
  
  for(i in 1:length(pDataFactors)){
    gg_color <- c(gg_color_hue(length(levels(pDataFactors[,i]))))[pDataFactors[,i]]
    matrix  <- plyr::mapvalues(pDataFactors[,i],
                               from = levels(pDataFactors[,i]),
                               to = unique(gg_color)) %>% 
      as.matrix()
    ColSideColors[,i] <- matrix 
    legendNames <- c(legendNames, "", levels(pDataFactors[,i]))
    legendColors <- c(legendColors, "white", unique(gg_color))
  }
  colnames(ColSideColors) <- names(pDataFactors)
  
  print(glue::glue("Obtaining smoothed methylation values..."))
  smoothed <- data.frame(getMeth(BSseq = bsseq, regions = regions, type = "smooth", what = "perRegion"))
  
  print(glue::glue("Tidying for heatmap of HCA..."))
  matrix <- as.matrix(smoothed)
  matrix <- matrix[,]*100
  
  # Subtract the mean methylation for each row/DMR (Trick to toggle)
  data <- sweep(matrix, 1, rowMeans(matrix))
  
  data <- as.matrix(data)
  colnames(data) <- groups
  
  print(glue::glue("Plotting heatmap of modified HCA Z-scores (% mCG/CG - mean)..."))
  source("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")
  heatmap.3(data,
            Rowv= as.dendrogram(hclust(dist(data))),
            scale = c("row"),
            Colv = TRUE,
            col = rev(brewer.pal(11, name = "RdBu")),
            margins = c(10,10),
            trace = "none",
            main = glue::glue("{length(regions)} Differentially Methylated Regions"),
            labRow = NA,
            #keysize = 0.85,
            #key.par = list(cex=0.5),
            KeyValueName = "Z-score(% mCG/CG - mean)",
            ColSideColors = ColSideColors,
            ColSideColorsSize = 2,
            ...
  )
  
  par(lend = 1)
  legend("topright",
         legend = legendNames,
         fill = legendColors,
         border = FALSE,
         bty = "n",
         y.intersp = 0.7,
         cex = 0.7)

}
