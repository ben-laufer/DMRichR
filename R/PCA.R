
#' smoothPCA
#' @description Provides individual smoothed methylation values for genomic ranges objects using bsseq
#' @param matrix Matrix of transposed individual methylation values
#' @param title Character string of title for plot and pdf
#' @return PCA plot
#' @require ggbiplot
#' @export PCA
PCA <- function(matrix = matrix,
                title = title){

  cat("\n[DMRichR] PCA \t\t\t\t\t\t", format(Sys.time(), "%d-%m-%Y %X"))
  message("Performing PCA...")
  data.pca <- prcomp(matrix, center = TRUE, scale. = TRUE)
  plot(data.pca, type = "l")
  print(summary(data.pca))

  message("Plotting PCA...")
  PCA <- ggbiplot(data.pca,
                  obs.scale = 1,
                  var.scale = 1,
                  groups = group,
                  ellipse = TRUE,
                  circle = FALSE,
                  var.axes = FALSE,
                  choices = 1:2) +
    scale_color_discrete(name = '') +
    theme_bw(base_size = 25) +
    geom_point(aes(colour = group), size = 4) +
    theme(legend.direction = 'vertical',
          legend.position = c(0.125, 0.1), # Change legend position
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 18),
          panel.grid.major = element_blank(),
          panel.border = element_rect(color = "black", size = 1.25),
          axis.ticks = element_line(size = 1.25),
          legend.key = element_blank(),
          panel.grid.minor = element_blank()) +
    guides(col=guide_legend(ncol=2)) +
    ggtitle(title) + # Change title
    theme(plot.title = element_text(hjust = 0.5))

  message("Saving PCA plot...")
  ggsave(paste(title,".pdf", sep = ""), plot = PCA, device = NULL)
  return(PCA)
}
