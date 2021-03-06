##' @title Plot the prior and posterior distribution of shifts
##'
##' @description Generates a barplot of the prior and posterior distributions
##'     of the number of shifts.
##'
##' @param mcmc A dataframe of the mcmc_out file from a \code{BAMM} run, or
##'     the filename.
##' @param expectedNumberOfShifts Expected number of shifts under the prior.
##' @param burnin The fraction of samples to discard as burn-in.
##' @param priorCol Color for the prior distribution.
##' @param postCol Color for the posterior distribution.
##' @param legendPos Placement of the legend, see \code{\link{legend}}.
##' @param \dots Additional parameters that are passed to
##'     \code{\link{barplot}}.
##'
##' @return Invisibly returns a matrix with the probability of each shift
##'     number under the prior and the posterior.
##'
##' @author Pascal Title
##'
##' @examples
##' data(mcmc.whales)
##' plotPrior(mcmc.whales, expectedNumberOfShifts = 1, burnin = 0.15)
##' @export
plotPrior <- function(mcmc, expectedNumberOfShifts = 1, burnin = 0.15, priorCol = 'light blue', postCol = 'red', legendPos = 'topright', ...) {
	
	if (!any(inherits(mcmc, c('character', 'data.frame', 'matrix')))) {
		stop('mcmc must be either a dataframe or the path to the mcmc_out file.')
	}
	
	if (is.character(mcmc)) {
		mcmc <- read.csv(mcmc, stringsAsFactors = FALSE)
	}
	
	#drop burnin
	mcmc2 <- mcmc[floor(burnin * nrow(mcmc)):nrow(mcmc),]
	
	#get prior distribution of shifts
	obsK <- seq(from = 0, to = max(mcmc2[,"N_shifts"]), by = 1)
	prior <- sapply(obsK, prob.k, poissonRatePrior = 1/expectedNumberOfShifts)
	prior <- data.frame(N_shifts = obsK, prob = prior)
	
	#get posterior distribution of shifts
	posterior <- sapply(obsK, function(x) length(which(mcmc2[,'N_shifts'] == x))) / nrow(mcmc2)
	names(posterior) <- obsK
	posterior <- data.frame(N_shifts = names(posterior), prob = posterior)

	barplot(prior[,2], names.arg = prior[,1], ylim = c(0, max(c(prior[,2], posterior[,2]))), border = 'black', col = priorCol, xlab = 'n shifts', ...)
	barplot(posterior[,2], add = TRUE, border = 'black', col = BAMMtools::transparentColor(postCol, 0.4), axes=FALSE)
	legend(x = legendPos, y = NULL, legend = c('prior','posterior'), fill = c(priorCol, BAMMtools::transparentColor(postCol, 0.4)), bty = 'n', cex=1.5)
	
	invisible(cbind(N_shifts = prior$N_shifts, priorProbs = prior$prob, postProbs = posterior$prob))	
}


prob.k <- function(k, poissonRatePrior) {
	Denom <- (poissonRatePrior + 1) ^ (k + 1)
	Prob <- poissonRatePrior / Denom
	return(Prob)
}
