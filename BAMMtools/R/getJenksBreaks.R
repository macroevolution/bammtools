###############################
#
#	getJenksBreaks <- function(...)
#	Function returns the Jenks natural breaks for a vector of values
#
#		var = numeric vector
#		k = number of categories
#		subset = if not NULL, this is the number of regularly sampled values to be taken from var


##' @title Jenks natural breaks classification
##'
##' @description Given a vector of numeric values and the number of desired
##'     breaks, calculate the optimum breakpoints using Jenks natural
##'     breaks optimization.
##'
##' @param var Numeric vector.
##' @param k Number of breaks.
##' @param subset Number of regularly spaced samples to subset from
##'     \code{var}. Intended to improve runtime for large datasets. If
##'     \code{NULL}, all values are used.
##'
##' @details \code{getJenksBreaks} is called by
##'     \code{\link{assignColorBreaks}}.
##' 
##'     The values in \code{var} are binned into \code{k+1} categories,
##'     according to the Jenks natural breaks classification method. This
##'     method is borrowed from the field of cartography, and seeks to
##'     minimize the variance within categories, while maximizing the variance
##'     between categories. If \code{subset = NULL}, all values of \code{var}
##'     are used for the optimization, however this can be a slow process with
##'     very large datasets. If \code{subset} is set to some number, then 
##'     \code{subset} regularly spaced values of \code{var} will be sampled.
##'     This is slightly less accurate than when using the entirety of
##'     \code{var} but is unlikely to make much of a difference. If 
##'     \code{subset} is defined but \code{length(var) < subset}, then
##'     \code{subset} has no effect.  
##' 
##'     The Jenks natural breaks method was ported to C from code found in the
##'     classInt R package.
##'
##' @return A numeric vector of intervals.
##'
##' @author Pascal Title
##'
##' @seealso See \code{\link{assignColorBreaks}} and
##'     \code{\link{plot.bammdata}}.
##'
##' @examples
##' # load whales dataset
##' data(whales, events.whales)
##' ed <- getEventData(whales, events.whales, burnin=0.25, nsamples=500)
##' 
##' # for demonstration purposes, extract the vector of speciation rates
##' ed <- dtRates(ed, tau=0.01)
##' vec <- ed$dtrates$rates[[1]]
##' 
##' # Return breaks for the binning of speciation rates into 65 groups
##' # yielding 64 breaks
##' getJenksBreaks(vec, 64)
##' @keywords graphics
##' @export
getJenksBreaks <- function(var, k, subset = NULL) {
	k <- k - 1;
	
	#if more breaks than unique values, segfault, so avoid
	if (k > length(unique(var))) {
		k <- length(unique(var));
	}
	brks <- rep(1, k + 1);
	
	#if requested, regularly sample subset values
	if (!is.null(subset)) {
		if (length(var) > subset) {
			ind <- c(seq(from=1, to=length(var), by=floor(length(var)/subset)), length(var));
			var <- var[ind];
		}
	}
	
	d <- sort(var);
	length_d <- length(d);
	return(.C("jenksBrks", as.double(d), as.integer(k), as.integer(length_d), as.double(brks))[[4]]);
}

