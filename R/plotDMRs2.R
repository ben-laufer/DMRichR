#' Plot Differentially Methylated Regions
#'
#' Generates trace plots of methylation proportions by genomic position. Modified
#' from the plotDMRs function in the dmrseq R package (\code{\link[dmrseq]{plotDMRs}}).
#' 
#' Creates aesthetially pleasing DMR plots. By default will plot individual
#' points with size proportional to coverage, along with a smoothed line 
#' for each sample. If BSseq object is smoothed, it will use these values. 
#' Otherwise, raw methylation values will be smoothed with loess. Elements 
#' will be colored by biological condition (\code{label}). Also has 
#' functionality to add annotations below the main plot (CpG category, genes) 
#' if \code{annoTrack} is specified.
#' 
#' @param BSseq An object of class BSseq. Can be smoothed or not smoothed.
#' 
#' @param regions A data.frame containing the DMRs (output from the main
#' \code{dmrseq} function).
#' 
#' @param testCovariate integer value or vector indicating which of columns of
#'  \code{pData(bs)} contains the covariate of interest. 
#'  This is used to construct the sample labels and colors (unless this is
#'  over-ridden by specifying \code{label}).
#'  
#' @param extend Describes how much the plotting region should be extended in
#'  either direction. The total width of the plot is equal to the width of the 
#'  region plus twice extend.
#' 
#' @param main The plot title. The default is to construct a title with 
#' information about which genomic region is being plotted.
#' 
#' @param addRegions A set of additional regions to be highlighted on the 
#' plots. Same format as the \code{regions} argument.
#' 
#' @param annoTrack a \code{GRangesList} object with two elements returned
#' by \code{\link{getAnnot}}. The first
#' contains CpG category information in the first element (optional)
#' coding gene sequence information in the second element (optional).
#' At least one of these elements needs to be non-null in order for 
#' any annotation to be plotted, but it is not necessary to contain
#' both.
#' 
#' @param col The color of the methylation estimates. It is recommended to 
#' leave this value as default (NULL), and specify a value of 
#' \code{testCovariate} to indicate which column of \code{pData(bs)}
#' to use as a factor for coloring the points and lines of the plot.
#' Alternatively, you can specify particular colors by 
#' passing this information through the \code{pData} slot of the
#' object \code{BSseq} (a \code{data.frame} that houses metadata). To do 
#' so, place the color value for each sample in a column titled \code{col}, 
#' and leave this argument as its default value of NULL. Alternatively,
#' you may specify a vector of color values (one for each sample), but 
#' you *must* make sure that this vector is in the same order as the samples
#' are in the BSseq object. If NULL and no \code{col} column is found in 
#' \code{pData}, then estimates are plotted in black for all samples.
#' 
#' @param lty The line type of the methylation estimates. It is recommended to
#' pass this information through the \code{pData} slot of the
#' object \code{BSseq} (a \code{data.frame} that houses metadata). To do 
#' so, place the line type value for each sample in a column titled \code{lty}, 
#' and leave this argument as its default value of NULL. Alternatively,
#' you may specify a vector of line type values (one for each sample), but 
#' you *must* make sure that this vector is in the same order as the samples
#' are in the BSseq object. If NULL and no \code{lty} column is found in 
#' \code{pData}, then estimates are plotted with \code{lty=1} for all samples.
#' 
#' @param lwd The line width of the methylation estimates. It is recommended to
#' pass this information through the \code{pData} slot of the
#' object \code{BSseq} (a \code{data.frame} that houses metadata). To do 
#' so, place the line width value for each sample in a column titled \code{lwd},
#' and leave this argument as its default value of NULL. Alternatively,
#' you may specify a vector of line width values (one for each sample), but 
#' you *must* make sure that this vector is in the same order as the samples
#' are in the BSseq object. If NULL and no \code{lwd} column is found in 
#' \code{pData}, then estimates are plotted with \code{lwd=1} for all samples.
#' 
#' @param label The condition/population labels for the plot legend. If NULL
#' (default) this is taken from the \code{testCovariate} column of 
#' \code{pData}. Alternatively, you can pass in labels by 
#' adding this information through the \code{pData} slot of the
#' object \code{BSseq} (a \code{data.frame} that houses metadata). To do 
#' so, place the labels for each sample in a column titled \code{label}, 
#' and leave this argument as its default value of NULL.
#' You may instead specify an arbitrary vector of labels (one for each sample),
#' but be aware that you *must* make sure that this vector is in the same order 
#' as the samples are in the BSseq object. If NULL, and \code{testCovariate} is
#' also NULL and no \code{label} column is found in 
#' \code{pData}, then no legend is created.
#' 
#' @param mainWithWidth logical value indicating whether the default title
#' should include information about width of the plot region.
#' 
#' @param regionCol The color used for highlighting the region.
#' 
#' @param addTicks logical value indicating whether tick marks showing the
#'  location of methylation loci should be added. Default is TRUE.
#' 
#' @param addPoints logical value indicating whether the individual 
#' methylation estimates be plotted as points.
#' 
#' @param pointsMinCov The minimum coverage a methylation loci need in 
#' order for the raw methylation estimates to be plotted. Useful for filtering
#' out low coverage loci. Only used if addPoints = TRUE. Default value is 1 
#' (no filtering).
#' 
#' @param highlightMain logical value indicating whether the plot region 
#' should be highlighted.
#' 
#' @param qval logical value indicating whether the region FDR estimate  
#' (q-value) should be displayed in the plot title. The value is extracted 
#' from the \code{regions} argument.
#' 
#' @param stat logical value indicating whether the region statistic 
#' should be displayed in the plot title. The value is extracted from the
#' \code{regions} argument.
#' 
#' @param verbose logical value indicating whether progress messages 
#' should be printed to the screen.
#' 
#' @param includeYlab a logical indicating whether to include the Y axis
#'  label 'Methylation' (useful to turn off if combining multiple region
#'  figures and you do not want to include redundant y axis label information)
#'  
#' @param compareTrack a named GenomicRangesList item that contains up to four
#' custom tracks (GenomicRanges objects) which will be plotted below the region.
#' Only one of `compareTrack` or `annoTrack` can be specified since there is
#' only for plotting either the built in GpG category and exon tracks, *or* a
#' custom set of tracks.
#' 
#' @param labelCols a character vector with names of the mcols slot of the 
#' GenomicRanges items in `compareTrack'. Only used if plotting custom
#' tracks using the `compareTrack' argument. If specified, the (first) value
#' in that column is printed along with a label that includes the name of the
#' list item. If NULL (default), just the name of the track is printed.
#' 
#' @param horizLegend logical indicating whether the legend should be 
#' horizontal instead of vertical (default FALSE). This is useful if you need
#' to plot many labels and want to preserve whitespace.
#' 
#' @param addLines logical indicating whether to plot smooth lines between 
#' points. Default is true. Can be useful to turn this off for very small 
#' regions.
#' 
#' @return None (generates a plot)
#' 
#' @import dmrseq
#' @import bsseq
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices hcl rainbow
#' @importFrom graphics arrows
#' 
#' @references \code{\link[dmrseq]{plotDMRs}}
#' 
#' @export plotDMRs2

plotDMRs2 <- function (BSseq, regions = NULL, testCovariate = NULL, extend = (end(regions) - start(regions) + 1)/2, 
                       main = "", addRegions = regions, annoTrack = NULL, col = NULL, lty = NULL, lwd = NULL, 
                       label = NULL, mainWithWidth = TRUE, regionCol = dmrseq:::.alpha("#C77CFF", 0.2), addTicks = TRUE, 
                       addPoints = TRUE, pointsMinCov = 1, highlightMain = FALSE, qval = TRUE, stat = TRUE, 
                       verbose = TRUE, includeYlab = TRUE, compareTrack = NULL, labelCols = NULL, 
                       horizLegend = FALSE, addLines = TRUE) {
        
        if (!addLines && !addPoints){ 
                stop("At least one of addLines or addPoints must be true")
        }
        if (verbose){
                message("[plotDMRs] Plotting ", length(regions$L), " DMRs")
        }
        if (addLines & verbose){
                if(bsseq::hasBeenSmoothed(BSseq)){
                        message("[plotDMRs] Smoothing BSmooth methylation values with loess for line plots")
                }
                else {
                        message("[plotDMRs] Smoothing raw methylation values with loess for line plots")
                }
        }
        if (!is.null(regions)) {
                if (is(regions, "data.frame")) {
                        gr <- data.frame2GRanges(regions, keepColumns = FALSE)
                }
                else {
                        gr <- regions
                }
                if (!is(gr, "GRanges")){ 
                        stop("'regions' needs to be either a 'data.frame' or a 'GRanges' ")
                }
        }
        else {
                gr <- granges(BSseq)
        }
        gr <- suppressWarnings(resize(gr, width = 2 * extend + width(gr), fix = "center"))
        overlap <- findOverlaps(BSseq, gr)
        overlapExtend <- NULL
        for(i in 1:length(unique(overlap@to))){
                temp <- overlap@from[overlap@to == unique(overlap@to)[i]]
                temp <- c(temp[1] - 1, temp, temp[length(temp)] + 1) # Add adjacent CpGs
                overlapExtend <- c(overlapExtend, temp)
        }
        overlapExtend <- overlapExtend[overlapExtend > 0 & overlapExtend <= length(BSseq)]
        BSseq <- BSseq[overlapExtend,]
        
        if (!is.null(annoTrack) && !is.null(compareTrack)){
                stop("Choose either annoTrack or compareTrack; can't plot both")
        }
        if (length(start(BSseq)) == 0){
                stop("No overlap between BSseq data and regions")
        }
        if (!is.null(main) && length(main) != length(gr)){
                main <- rep(main, length = length(gr))
        }
        if (length(extend) == 1) {
                extend <- rep(extend, length(gr))
        }
        if (!is.null(testCovariate)) {
                coeff <- seq(2, (1 + length(testCovariate)))
                testCov <- as.character(pData(BSseq)[, testCovariate])
                if (length(unique(testCov)) > 2 && !is.numeric(testCov) && length(coeff) == 1) {
                        coeff <- seq(coeff, coeff + length(unique(as.character(testCov))) - 2)
                }
                design <- model.matrix(~testCov)
                if (is.null(col) && !("col" %in% names(pData(BSseq)))) {
                        cov.unique <- unique(design[, coeff, drop = FALSE])
                        ncol <- nrow(cov.unique)
                        colors <- gg_color_hue(ncol)
                        if (ncol == 2) {
                                colors <- c("mediumblue", "deeppink1")
                        }
                        colors <- cbind(cov.unique, colors[rank(as.numeric(rowSums(cov.unique)), ties = "first")])
                        colmat <- colors[, -ncol(colors), drop = FALSE]
                        colmat <- apply(colmat, 2, as.numeric)
                        z <- colors[, ncol(colors)][match(data.frame(t(design[,coeff, drop = FALSE])), data.frame(t(colmat)))]
                        pData(BSseq)$col <- as.character(z)
                }
                if (is.null(label) && !("label" %in% names(pData(BSseq)))) {
                        pData(BSseq)$label <- paste0(pData(BSseq)[, testCovariate])
                }
        }
        if (!is.null(label) || "label" %in% names(pData(BSseq))) {
                if (!is.null(label)) {
                        labs <- label
                }
                else {
                        labs <- pData(BSseq)[["label"]]
                }
                if (horizLegend) {
                        wiggle <- max(nchar(labs)) * 0.5
                }
                else {
                        wiggle <- length(unique(labs)) * 0.9
                }
                opar <- par(mar = c(0, 4.1, 0.5, wiggle), oma = c(0, 0.5, 1, 1.5))
        }
        else {
                opar <- par(mar = c(0, 4.1, 0.5, 0), oma = c(0, 0.5, 1, 1.5))
        }
        on.exit(par(opar))
        for (ii in seq(along = gr)) {
                if (verbose && ii%%100 == 0) {
                        cat(sprintf("..... Plotting region %d (out of %d)\n", 
                                    ii, nrow(regions)))
                }
                .plotSingleDMR2(BSseq = BSseq, region = regions[ii, ], 
                                extend = extend[ii], main = main[ii], col = col, 
                                lty = lty, lwd = lwd, label = label, addRegions = addRegions, 
                                regionCol = regionCol, mainWithWidth = mainWithWidth, 
                                annoTrack = annoTrack, addTicks = addTicks, addPoints = addPoints, 
                                pointsMinCov = pointsMinCov, highlightMain = highlightMain, 
                                qval = qval, stat = stat, includeYlab = includeYlab, 
                                compareTrack = compareTrack, labelCols = labelCols, 
                                horizLegend = horizLegend, addLines = addLines)
        }
}

.plotSingleDMR2 <- function (BSseq, region = NULL, extend = 0, main = "", addRegions = NULL, annoTrack = NULL, 
                             col = NULL, lty = NULL, lwd = NULL, label = NULL, mainWithWidth = TRUE, 
                             regionCol = dmrseq:::.alpha("orchid1", 0.2), addTicks = TRUE, addPoints = FALSE, 
                             pointsMinCov = 5, highlightMain = FALSE, qval = NULL, stat = NULL, 
                             includeYlab = TRUE, compareTrack = NULL, labelCols = NULL, horizLegend = FALSE, 
                             addLines = TRUE){
        # Modified from the dmrseq R package
        if (!is.null(annoTrack) || !is.null(compareTrack)) {
                layout(matrix(seq_len(2), ncol = 1), heights = c(2, 1.25))
        }
        else {
                layout(matrix(seq_len(2), ncol = 1), heights = c(2, 0.2))
        }
        .dmrPlotSmoothData2(BSseq = BSseq, region = region, extend = extend, 
                            addRegions = addRegions, col = col, lty = lty, lwd = lwd, 
                            label = label, regionCol = regionCol, addTicks = addTicks, 
                            addPoints = addPoints, pointsMinCov = pointsMinCov, highlightMain = highlightMain, 
                            includeYlab = includeYlab, horizLegend = horizLegend, 
                            addLines = addLines)
        gr <- dmrseq:::bsseq.bsGetGr(BSseq, region, extend)
        if (!is.null(main)) {
                if (qval && stat) {
                        qval <- round(region$qval, 4)
                        stat <- round(region$stat, 3)
                        main <- dmrseq:::.dmrPlotTitle(gr = region, extend = extend, 
                                              main = main, mainWithWidth = mainWithWidth, qval = qval, 
                                              stat = stat)
                }
                else if (stat) {
                        stat <- round(region$stat, 3)
                        main <- dmrseq:::.dmrPlotTitle(gr = region, extend = extend, 
                                              main = main, mainWithWidth = mainWithWidth, stat = stat)
                }
                else if (qval) {
                        qval <- round(region$qval, 4)
                        main <- dmrseq:::.dmrPlotTitle(gr = region, extend = extend, 
                                              main = main, mainWithWidth = mainWithWidth, qval = qval)
                }
                else {
                        main <- dmrseq:::.dmrPlotTitle(gr = region, extend = extend, 
                                              main = main, mainWithWidth = mainWithWidth)
                }
                mtext(side = 3, text = main, outer = FALSE, cex = 1, line = 0, padj = 0)
        }
        if (!is.null(annoTrack)) {
                .dmrPlotAnnotations2(gr, annoTrack)
        }
        else if (!is.null(compareTrack)) {
                dmrPlotComparisons(gr, compareTrack, labelCols = labelCols)
        }
}

.dmrPlotAnnotations2 <- function (gr, annoTrack){
        # Modified from the dmrseq R package
        if (!is(annoTrack, "SimpleGenomicRangesList")) 
                stop("'annoTrack' needs to be a 'SimpleGenomicRangesList'")
        plot(start(gr), 1, type = "n", xaxt = "n", yaxt = "n", bty = "n", 
             ylim = c(0.25, length(annoTrack) - 0.1), xlim = c(start(gr), end(gr)), xlab = "", ylab = "")
        vars <- list(Island = "Island   ", Shore = "Shore   ", Shelf = "Shelf   ", 
                     OpenSea1 = "Open ", OpenSea2 = "Sea")
        cols <- c("forestgreen", "goldenrod2", "dodgerblue", "blue3", "blue3")
        for (i in seq_len(length(vars))) {
                tmpvars <- vars
                tmpvars[-i] <- paste("phantom('", tmpvars[-i], "')", sep = "")
                expr <- paste(tmpvars, collapse = "*")
                text(start(gr), 1.5, parse(text = expr), col = cols[i], adj = c(0, 1), cex = 1)
        }
        lapply(seq(along = annoTrack), function(ii) {
                jj <- length(annoTrack) + 1 - ii
                ir <- subsetByOverlaps(annoTrack[[ii]], gr)
                start(ir) <- pmax(start(ir), start(gr))
                end(ir) <- pmin(end(ir), end(gr))
                if (ii == 2) {
                        jj <- jj - 0.15
                }
                if (length(ir) > 0) {
                        if (ii == 2) {
                                colourCount <- length(unique(ir$symbol))
                                getPalette <- colorRampPalette(RColorBrewer::brewer.pal(max(length(unique(ir$symbol)), 3), "Dark2"))
                                color.pal <- getPalette(colourCount)
                                names(color.pal) <- unique(ir$symbol)
                                map <- match(ir$symbol, names(color.pal))
                                color <- color.pal[map]
                                bord <- "black"
                        }
                        else if (ii == 1) {
                                color <- ir$type
                                colvec <- c("blue3", "dodgerblue", "goldenrod2", "forestgreen")
                                names(colvec) <- c("inter", "shelves", "shores", "islands")
                                
                                for (ucol in unique(color)){
                                        ix <- agrep(ucol, names(colvec))
                                        if (length(ix) > 1)
                                                stop("Ambiguous CpG annotation labels")
                                        
                                        color[color == ucol] <- colvec[ix]
                                }
                                bord <- color
                                rect(start(ir), jj - 0.4, end(ir), jj - 0.15, col = color, border = bord)
                        }
                        if (ii == 2) {
                                lastPos <- rep(NA, length(unique(ir$symbol)) - 1)
                                used <- NULL
                                for (k in seq_len(length(unique(ir$symbol)))) {
                                        irk <- ir[ir$symbol == unique(ir$symbol)[k], ]
                                        rect(min(start(irk)), jj - 0.08, max(end(irk)), 
                                             jj + 0.08, col = dmrseq:::.alpha("black", 0.1), border = dmrseq:::.alpha("black", 0.1))
                                        if (!(unique(ir$symbol)[k] %in% used)) {
                                                rwidth <- end(gr) - start(gr)
                                                sg <- pmin(start(irk), end(irk))
                                                eg <- pmax(start(irk), end(irk))
                                                gwidth <- min(max(eg), end(gr)) - max(min(sg), start(gr))
                                                textPos <- max(min(sg), start(gr)) + gwidth/2
                                                jj.orig <- jj
                                                if (sum(!is.na(lastPos)) > 0) {
                                                        separation <- (textPos - lastPos[k - 1])/rwidth
                                                        if (abs(separation) <= 0.2 && k < 3) {
                                                                jj <- jj - 0.16
                                                        }
                                                        else {
                                                                separation <- min(abs((textPos - lastPos)/rwidth), na.rm = TRUE)
                                                                if (abs(separation) <= 0.2) {
                                                                        jj <- jj - 0.16
                                                                }
                                                        }
                                                }
                                                lastPos[k] <- textPos
                                                text(textPos, jj - 0.375, labels = unique(irk$symbol), 
                                                     cex = 1, col = unique(color)[k], font = 3)
                                                jj <- jj.orig
                                                used <- c(used, unique(ir$symbol)[k])
                                        }
                                        rect(sg, jj - 0.22, eg, jj + 0.22, col = unique(color)[k], border = bord)
                                }
                        }
                }
                mtext(names(annoTrack)[ii], side = 2, at = jj * 0.75 + 0.22, las = 1, line = 1)
        })
}

.dmrPlotSmoothData2 <- function (BSseq, region, extend, addRegions, col, lty, lwd, label, regionCol, addTicks, 
                                 addPoints, pointsMinCov, highlightMain, includeYlab = TRUE, horizLegend, 
                                 addLines = TRUE){
        # Modified from the dmrseq R package
        gr <- dmrseq:::bsseq.bsGetGr(BSseq, region, extend)
        overlap <- findOverlaps(BSseq, gr)
        overlapExtend <- overlap@from[overlap@to == unique(overlap@to)[1]] # Only 1 region
        overlapExtend <- c(overlapExtend[1] - 1, overlapExtend, overlapExtend[length(overlapExtend)] + 1) # Add adjacent CpGs
        BSseq <- BSseq[overlapExtend,]
        sampleNames <- sampleNames(BSseq)
        names(sampleNames) <- sampleNames
        positions <- start(BSseq)
        rawPs <- as.matrix(bsseq::getMeth(BSseq, type = "raw"))
        if(bsseq::hasBeenSmoothed(BSseq)){
                smoothPs <- as.matrix(bsseq::getMeth(BSseq, type = "smooth"))
        }
        coverage <- as.matrix(bsseq::getCoverage(BSseq))
        colEtc <- dmrseq:::.dmrGetMeta(object = BSseq, col = col, lty = lty, lwd = lwd, label = label)
        if (includeYlab) {
                yl <- "Methylation"
        }
        else {
                yl <- ""
        }
        plot(positions[1], 0.5, type = "n", xaxt = "n", yaxt = "n", ylim = c(-0.06, 1), 
             xlim = c(start(gr), end(gr)), xaxs = "i", xlab = "", ylab = yl)
        axis(side = 2, at = c(0, 0.25, 0.5, 0.75, 1), labels = c(0, 25, 50, 75, 100), las = 1)
        if (addTicks){ 
                suppressWarnings(rug(positions, ticksize = 0.06))
        }
        if (is.list(addRegions) && !is.data.frame(addRegions)) {
                if (length(addRegions) > 2) {
                        stop("Only two sets of regions can be highlighted")
                }
                if (length(regionCol) == 1) {
                        regionCol <- c(regionCol, dmrseq:::.alpha("blue", 0.2))
                }
                dmrseq:::bsseq.bsHighlightRegions(regions = addRegions[[1]], gr = gr, ylim = c(0, 1), 
                                         regionCol = regionCol[1], highlightMain = highlightMain)
                dmrseq:::bsseq.bsHighlightRegions(regions = addRegions[[2]], gr = gr, ylim = c(0, 1), 
                                         regionCol = regionCol[2], highlightMain = highlightMain)
        }
        else {
                dmrseq:::bsseq.bsHighlightRegions(regions = addRegions, gr = gr, ylim = c(0, 1), 
                                         regionCol = regionCol, highlightMain = highlightMain)
        }
        if (addPoints) {
                for (sampIdx in seq_len(ncol(BSseq))) {
                        dmrseq:::.dmrPlotPoints(x = positions[2:(length(positions) - 1)], 
                                                y = rawPs[2:(nrow(rawPs) - 1), sampIdx], 
                                                z = coverage[2:(nrow(coverage) - 1), sampIdx], 
                                                col = colEtc$col[sampIdx], pointsMinCov = pointsMinCov, 
                                                maxCov = quantile(coverage, 0.95))
                }
        }
        if (addLines) {
                for (sampIdx in seq_len(ncol(BSseq))) {
                        if (sum(!is.na(rawPs[, sampIdx])) > 1) {
                                if(bsseq::hasBeenSmoothed(BSseq)){
                                        .dmrPlotSmoothLines(x = positions, y = smoothPs[,sampIdx], 
                                                            z = rep(1, sum(!is.na(smoothPs[,sampIdx]))),
                                                            col = dmrseq:::.makeTransparent(colEtc$col[sampIdx], 200), 
                                                            lwd = colEtc$lwd[sampIdx], lty = colEtc$lty[sampIdx],
                                                            spn = max(1 - sum(!is.na(smoothPs[,sampIdx])) / 50, 0.1))
                                }
                                else {
                                        .dmrPlotSmoothLines(x = positions, y = rawPs[,sampIdx], z = coverage[,sampIdx],
                                                            col = dmrseq:::.makeTransparent(colEtc$col[sampIdx], 200), 
                                                            lwd = colEtc$lwd[sampIdx], lty = colEtc$lty[sampIdx],
                                                            spn = max(1 - sum(!is.na(rawPs[,sampIdx])) / 160, 0.75)) 
                                }
                        }
                }
        }
        if (sum(!is.na(colEtc$label)) == length(colEtc$label)) {
                .dmrPlotLegend2(plotRange = c(start(gr), end(gr)), colEtc$col, colEtc$label, horizLegend)
        }
        box() # Plot the border in front
}

.dmrPlotSmoothLines <- function (x, y, z, col, lwd, lty, spn){
        # Modified from the dmrseq R package
        y[y == 1] <- 0.99
        y[y == 0] <- 0.01
        logit <- function(p) {
                log(p/(1 - p))
        }
        inv.logit <- function(l) {
                exp(l)/(1 + exp(l))
        }
        if (length(x) >= 10) {
                loess_fit <- suppressWarnings(loess(logit(y) ~ x, weights = z, span = spn))
                xl <- seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE), length.out = 500)
                lines(xl, inv.logit(predict(loess_fit, xl)), col = col, lwd = lwd, lty = lty)
        }
        else {
                lines(x, y, col = col, lwd = lwd, lty = lty)
        }
}

.dmrPlotLegend2 <- function (plotRange, col, label, horizLegend){
        # Modified from the dmrseq R package
        numUnique <- length(unique(paste0(col, label, sep = "")))
        if (numUnique < length(col)) {
                col <- unique(col)
                label <- unique(label)
        }
        if (!horizLegend) {
                for (lg in seq_len(length(label))) {
                        mtext(label[lg], side = 4, line = lg - 1, col = col[lg], cex = 0.9, las = 0)
                }
        }
        else {
                for (lg in seq_len(length(label))) {
                        mtext(label[lg], side = 4, line = 0.5, col = col[lg], cex = 0.9, las = 1, 
                              at = 1.02 - 0.08 * (lg - 1))
                }
        }
}
