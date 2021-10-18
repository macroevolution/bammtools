
# ratesHistogram <- function(...)
		# phylorates = a plot.bammdata object
		# plotBrks = boolean, should breaks be plotted
		# xlab = x-axis label
		# ylab = y-axis label
		# lwd = lwd for breaks lines
		# lty = lty for breaks lines
		# brksCol = color for breaks lines
		# ... additional arguments passed on to mtext for axis labels

##' @title Histogram of \code{BAMM} rate frequencies
##'
##' @description Plots a histogram of the frequency of rate values across the
##'     phylogeny.
##'
##' @param phylorates A saved \code{\link{plot.bammdata}} object.
##' @param plotBrks Boolean, should breaks be plotted over the histogram.
##' @param xlab x-axis label.
##' @param ylab y-axis label.
##' @param lwd Line width for breaks.
##' @param lty Line style for breaks.
##' @param brksCol Color of breaks lines.
##' @param \dots Additional arguments passed on to \code{mtext} for axis
##'     labels.
##'
##' @details With this function, a histogram is plotted that shows the
##'     frequency of rates present in the dataset. The color scheme plotted
##'     is taken from the saved \code{plot.bammdata} object that is the main
##'     input. Therefore, the mapping of colors to rates in the histogram
##'     corresponds exactly to what is plotted in the phylorate plot. If
##'     \code{plotBrks = TRUE}, then the color breaks used for the phylorates
##'     plot are shown. 
##'
##'     This function can be a useful tool for exploring different
##'     \code{plot.bammdata} options. Please see
##'     \url{http://bamm-project.org/colorbreaks.html} on the bamm-project
##'     website for more information on the utility of this function.
##'
##' @author Pascal Title
##'
##' @seealso \code{\link{plot.bammdata}}
##'
##' @references \url{http://bamm-project.org/}
##'
##' @examples
##' data(primates, events.primates)
##' ed <- getEventData(primates, events.primates, burnin=0.25, nsamples=500,
##'                    type = 'trait')
##'
##' # create phylorate plot with the jenks breaks method to generate output
##' phylorates <- plot(ed, breaksmethod='jenks', show = FALSE)
##'
##' ratesHistogram(phylorates, plotBrks = TRUE, xlab = 'trait rates')
##' @export
ratesHistogram <- function(phylorates, plotBrks = TRUE, xlab = 'speciation rate', ylab = 'density', lwd = 0.2, lty = 1, brksCol = 'black', ...) {
	
	if (!identical(names(phylorates), c("coords", "colorbreaks", "palette", "colordens"))) {
		stop("phylorates must be a saved plot.bammdata object.")
	}

	plot.new();
	x <- phylorates$colordens[,1];
	y <- phylorates$colordens[,2];
	plot.window(xlim = c(min(x), max(x)), ylim = c(0, max(y)));
	segments(x, y, x, 0, lend = 2, col = phylorates$colordens[,3], lwd = 3);
	axis(1, signif(seq(min(0, min(x)), max(x), length.out = 5), 2), xaxs = "i", cex.axis = 0.75, tcl = NA, mgp = c(0, 0.25, 0));
	axis(2, round(c(-1, seq(0, max(y), length.out = 4))), las = 1, yaxs = "i", cex.axis = 0.75, tcl = NA, mgp = c(0, 0.35, 0));
	
	mtext(xlab, side = 1, line = 1.5, ...);
	mtext(ylab, side = 2, line = 1.5, ...);
	
	#add breaks as vertical lines
	if (plotBrks) {
		abline(v = phylorates$colorbreaks, lwd = lwd, lty = lty, col = brksCol);
	}
}





