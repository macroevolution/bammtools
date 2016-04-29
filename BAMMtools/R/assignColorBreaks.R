##' @title Map macroevolutionary rates to colors
##'
##' @description Maps macroevolutionary rates to a set of \code{NCOLORS}.
##'
##' @param rates A numeric vector of phenotypic rates or a list of numeric
##'     vectors of speciation and extinction rates.
##' @param NCOLORS An integer number of colors to use for the mapping. Larger
##'     numbers do not necessarily result in smoother looking color ramps. The
##'     default is 64 and is probably sufficient for most purposes.
##' @param spex A character string. "s" means that speciation rates are used
##'     to make the map, "e" means that extinction rates are used. "netdiv"
##'     means that diversification rates are used. Ignored for \code{BAMM}
##'     trait data.
##' @param logcolor Logical. Should the natural logarithm of rates be used for
##'     the color map.
##' @param method Determines how the color breaks are created. See Details.
##' @param JenksSubset Number of regularly spaced samples to subset from
##'     \code{rates}. Only relevant when \code{method = "jenks"}. See Details.
##'
##' @details If \code{method = "quantile"} macroevolutionary rates are binned
##'     into \code{NCOLORS+1} percentiles and rates in each bin are mapped to
##'     a color determined by the \code{pal} argument in \code{plot.bammdata}.
##'     Alternatively, if \code{method = "linear"} macroevolutionary rates are
##'     binned into \code{NCOLORS+1} equal length intervals between the
##'     minimum and maximum. 
##'
##'     If \code{method = "jenks"}, macroevolutionary rates are binned into
##'     \code{NCOLORS+1} categories, according to the Jenks natural breaks
##'     classification method. This method is borrowed from the field of
##'     cartography, and seeks to minimize the variance within categories,
##'     while maximizing the variance between categories. 
##'
##'     The Jenks natural breaks method was ported to C from code found in the classInt R package. 
##'
##' @return A numeric vector of rate percentiles/intervals.
##'
##' @author Mike Grundler, Pascal Title
##'
##' @seealso \code{\link{plot.bammdata}}
##'
##' @examples
##' data(whales, events.whales)
##' ed <- getEventData(whales, events.whales, burnin = 0.2, nsamples = 500)
##'
##' ed <- dtRates(ed, 0.01)
##' colors <- assignColorBreaks(ed$dtrates$rates, spex="s") #speciation rates
##' #colors <- assignColorBreaks(ed$dtrates$rates[[1]]) 
##' #this also works for speciation rates
##'
##' plot(ed, colorbreaks = colors, spex="s")
##' colors <- assignColorBreaks(ed$dtrates$rates, spex="netdiv") 
##' #diversification rates
##'
##' #colors <- assignColorBreaks(ed$dtrates$rates[[1]] - ed$dtrates$rates[[2]]) 
##' #this also works for diversification rates
##'
##' plot(ed, colorbreaks = colors, spex="netdiv")
##' @keywords graphics
##' @export
assignColorBreaks <- function(rates, NCOLORS = 64, spex = "s", logcolor = FALSE, method = c("linear","quantile","jenks"), JenksSubset = NULL) {
	method = match.arg(method, c("linear", "quantile", "jenks"));
		if (mode(rates) == "numeric") {
		if (logcolor == FALSE) {
			if (method == "quantile") {
				bks <- quantile(rates, seq(0,1, length.out=(NCOLORS+1)));
			}
			if (method == 'jenks') {
				bks <- getJenksBreaks(rates, k=(NCOLORS + 1), subset = JenksSubset);
			}
			if (method == 'linear') {
				bks <- seq(min(rates), max(rates), length.out = (NCOLORS+1));
			}
		}
		else {
			if (method == "quantile") {
				bks <- quantile(log(rates), seq(0,1, length.out=(NCOLORS+1)));
			}
			if (method == 'jenks') {
				bks <- getJenksBreaks(log(rates), k=(NCOLORS + 1), subset = JenksSubset);
			}
			if (method == 'linear') {
				bks <- seq(min(log(rates)), max(log(rates)), length.out = (NCOLORS+1));
			}
		}	
	}
	else if (mode(rates) == "list") {
		if (tolower(spex) == "s") {
			if (logcolor == FALSE) {
				if (method == "quantile") {
					bks <- quantile(rates[[1]], seq(0,1, length.out=(NCOLORS+1)));
				}
				if (method == 'jenks') {
					bks <- getJenksBreaks(rates[[1]], k=(NCOLORS + 1), subset = JenksSubset);
				}
				if (method == 'linear') {
					bks <- seq(min(rates[[1]]), max(rates[[1]]), length.out = (NCOLORS+1));
				}
			}
			else {
				if (method == "quantile") {	
					bks <- quantile(log(rates[[1]]), seq(0,1, length.out=(NCOLORS+1)));
				}
				if (method == 'jenks') {
					bks <- getJenksBreaks(log(rates[[1]]), k=(NCOLORS + 1), subset = JenksSubset);
				}
				if (method == 'linear') {
					bks <- seq(min(log(rates[[1]])), max(log(rates[[1]])), length.out = (NCOLORS+1));
				}
			}
		}
		else if (tolower(spex) == "e") {
			if (logcolor == FALSE) {
				if (method == "quantile") {
					bks <- quantile(rates[[2]], seq(0,1, length.out=(NCOLORS+1)));
				}
				if (method == 'jenks') {
					bks <- getJenksBreaks(rates[[2]], k=(NCOLORS + 1), subset = JenksSubset);
				}
				if (method == 'linear') {
					bks <- seq(min(rates[[2]]), max(rates[[2]]), length.out = (NCOLORS+1));
				}
			}
			else {
				if (method == "quantile") {	
					bks <- quantile(log(rates[[2]]), seq(0,1, length.out=(NCOLORS+1)));
				}
				if (method == 'jenks') {
					bks <- getJenksBreaks(log(rates[[2]]), k=(NCOLORS + 1), subset = JenksSubset);
				}
				if (method == 'linear') {
					bks <- seq(min(log(rates[[1]])), max(log(rates[[1]])), length.out = (NCOLORS+1));
				}
			}
		}
		else if (tolower(spex) == "netdiv") {
			if (logcolor == FALSE) {
				if (method == "quantile") {	
					bks <- quantile(rates[[1]] - rates[[2]], seq(0,1, length.out=(NCOLORS+1)));
				}
				if (method == 'jenks') {
					bks <- getJenksBreaks(rates[[1]] - rates[[2]], k=(NCOLORS + 1), subset = JenksSubset);
				}
				if (method == 'linear') {
					bks <- seq(min(rates[[1]] - rates[[2]]), max(rates[[1]] - rates[[2]]), length.out = (NCOLORS+1));
				}
			}
			else { 
				z <- safeLog(rates[[1]] - rates[[2]]);
				if (method == "quantile") {
					bks <- quantile(z, seq(0,1, length.out=(NCOLORS+1)));
					#bks <- quantile(log(rates[[1]] - rates[[2]]), seq(0,1, length.out=(NCOLORS+1)))
				}
				if (method == 'jenks') {
					bks <- getJenksBreaks(z, k=(NCOLORS + 1), subset = JenksSubset);
				}
				if (method == 'linear') {
					bks <- seq(min(z), max(z), length.out = (NCOLORS+1));
					#bks <- seq(min(log(rates[[1]] - rates[[2]])), max(min(log(rates[[1]] - rates[[2]]))), length.out=(NCOLORS+1) );
				}
				attr(bks, "increment") <- attr(z, "increment");
				return (safeLog(bks, inverse = TRUE));
			}
		}
	}
	
	if (logcolor)
		return (exp(bks));
	return (bks)
}


safeLog <- function(x, inverse = FALSE) {
	
	if (inverse)
		y <- exp(x) - attr(x, "increment")
	else {
		y <- log(x + abs(min(x)) + 0.0001);
		attr(y, "increment") <- abs(min(x)) + 0.0001; 
	}	
	return (y);
}
