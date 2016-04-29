#############################################################
#
#	getTipRates(....)
#
# Returns a list with:
# lambda = matrix of tip rates where rows are species and columns are posterior samples, 
# mu if ephy$type == 'diversification',
# beta if ephy$type='trait',
# lambda.avg, mu.avg, beta.avg: named vector of average tip rates. 
# If returnNetDiv = TRUE, then a matrix and average vector for net div rates is returned.

##' @title Compute tip-specific macroevolutionary rates from \code{bammdata}
##'     object
##'
##' @description Return speciation, extinction, net diversification, or
##'     Brownian motion trait rates for all species in the phylogeny from
##'     \code{BAMM} output.
##'
##' @param ephy An object of class \code{bammdata}.
##' @param returnNetDiv Logical. If \code{TRUE}, then net diversification
##'     rates are returned, if \code{FALSE}, then both speciation and
##'     extinction rates are returned. If \code{ephy} is of type \code{trait},
##'     then this is ignored.
##' @param statistic Determines how the average tip rates should be
##'     calculated. Can be either \code{mean} or \code{median}.
##'
##' @return Returns a list with the following elements:
##'
##'     If \code{ephy} type is 'diversification':
##'     \itemize{
##'         \item{lambda} {A matrix of tip speciation rates with species as
##'             rows, and posterior samples as columns.}
##'         \item{mu} {A matrix of tip extinction rates with species as rows,
##'             and posterior samples as columns.}
##'         \item{lambda.avg} {A vector of average tip speciation rates,
##'             averaged with mean or median, depending on selected option for
##'             \code{statistic}. The vector is named with species names.}
##'         \item{mu.avg} {A vector of average tip extinction rates, averaged
##'             with mean or median, depending on selected option for
##'             \code{statistic}. The vector is named with species names.}
##'     }
##'
##'     If \code{ephy} type is 'diversification' and
##'     \code{returnNetDiv = TRUE}:
##'     \itemize{
##'         \item{netdiv} {A matrix of tip net diversification rates with
##'             species as rows, and posterior samples as columns.}
##'
##'         \item{netdiv.avg} {A vector of average tip net diversification
##'             rates, averaged with mean or median, depending on selected
##'             option for \code{statistic}. The vector is named with species
##'             names.}
##'     }
##'
##'     If \code{ephy} type is 'trait':
##'     \itemize{
##'         \item{beta} {A matrix of tip phenotypic rates with species as
##'             rows, and posterior samples as columns.}
##'         \item{beta.avg} {A vector of average tip phenotypic rates,
##'             averaged with mean or median, depending on selected option for
##'             \code{statistic}. The vector is named with species names.}
##'     }
##'
##' @author Pascal Title
##'
##' @seealso Requires an object of class \code{bammdata} as obtained with
##'     \code{\link{getEventData}}.
##'
##' @examples
##' data(whales, events.whales)
##' ephy <- getEventData(whales, events.whales, burnin=0.25, nsamples = 500)
##' 
##' # return a vector of average species-specific speciation rates.
##' meanlam <- getTipRates(ephy, returnNetDiv = FALSE,
##'                        statistic = 'mean')$lambda.avg
##' meanlam
##' 
##' # return a vector of median species-specific net diversification rates.
##' ndr <- getTipRates(ephy, returnNetDiv = TRUE,
##'                    statistic = 'median')$netdiv.avg
##' 
##' # Return mean species-specific speciation rates from all posterior 
##' # samples in the \code{bamm-data} object.
##' lam <- getTipRates(ephy, returnNetDiv = FALSE, statistic = 'mean')$lambda
##' rowMeans(lam)
##' @keywords models
##' @export
getTipRates <- function(ephy, returnNetDiv = FALSE, statistic = 'mean') {
	
	if (!'bammdata' %in% class(ephy)) {
		stop("Object ephy must be of class bammdata\n");
	}
	if (!statistic %in% c('mean','median')) {
		stop("statistic must be either 'mean' or 'median'.");
	}

	obj <- list();
	if (ephy$type == 'diversification') {
		if (returnNetDiv) {
			obj$netdiv <- do.call(cbind, ephy$tipLambda) - do.call(cbind, ephy$tipMu);
			rownames(obj$netdiv) <- as.phylo.bammdata(ephy)$tip.label;
			
			if (statistic == 'mean') {
				obj$netdiv.avg <- rowMeans(obj$netdiv);
			}
			if (statistic == 'median') {
				obj$netdiv.avg <- apply(obj$netdiv, 1, median);
			}
		}
		
		if (!returnNetDiv) {
			obj$lambda <- do.call(cbind, ephy$tipLambda);
			rownames(obj$lambda) <- as.phylo.bammdata(ephy)$tip.label;
		
			obj$mu <- do.call(cbind, ephy$tipMu);
			rownames(obj$mu) <- as.phylo.bammdata(ephy)$tip.label;
		
			if (statistic == 'mean') {
				obj$lambda.avg <- rowMeans(obj$lambda);
				obj$mu.avg <- rowMeans(obj$mu);
			}
			if (statistic == 'median') {
				obj$lambda.avg <- apply(obj$lambda, 1, median);
				obj$mu.avg <- apply(obj$mu, 1, median);
			}
		}
	}
	
	if (ephy$type == 'trait') {
		obj$beta <- do.call(cbind, ephy$tipLambda);
		rownames(obj$beta) <- as.phylo.bammdata(ephy)$tip.label;
		
		if (statistic == 'mean') {
			obj$beta.avg <- rowMeans(obj$beta);
		}
		if (statistic == 'median') {
			obj$betabeta.avg <- apply(obj$beta, 1, median);
		}
	}
	return(obj);
}









