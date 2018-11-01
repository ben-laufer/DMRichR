#' smoothHeatmap
#' @description Plot a heatmap of individual smoothed methylation value z scores for selected regions (i.e. significant DMRs)
#' @param bsseq Smoothed bsseq object
#' @param regions Genomic ranges object
#' @param names Ordered sample names
#' @param groups Ordered test covariate information for each sample
#' @param out Name of the text file to save in quotations
#' @return Saves a pdf image of the heatmap
#' @import bsseq
#' @import tidyverse
#' @import gplots
#' @export smoothHeatmap
smoothHeatmap <- function(regions = sigRegions,
                          bsseq = bs.filtered.bsseq,
                          groups = as.tibble(pData(bs.filtered.bsseq)) %>% pull(!!testCovariate),
                          out = "sig_individual_smoothed_DMR_methylation.txt",
                          ...){
  cat("\n[DMRichR] DMR heatmap \t\t\t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
  message("Smoothing...")
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
  pdf("heatmap.pdf", height = 8.5, width = 11)
  heatmap.2(data,
            Rowv= as.dendrogram(hclust(dist(data))),
            scale = c("row"),
            Colv = TRUE,
            col = rev(brewer.pal(11, name = "RdBu")),
            margins = c(10,10),
            trace = "none",
            main = paste(nrow(sig_indiv_smoothed),"Differentially Methylated Regions", sep = " "),
            labRow = NA,
            srtCol = 60,
            keysize = 0.85,
            key.par = list(cex=0.5),
            key.xlab = "Z-score(% mCG/CG - mean)",
            key.ylab = "Frequency",
            key.title = "",
            ...
  )
  dev.off()
}
