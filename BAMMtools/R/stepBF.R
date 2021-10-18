##' @title Identify the optimal number of shifts using Bayes factors
##'
##' @description stepBF is a function to determine the overall best fitting number of shifts 
##' using Bayes factor evidence.
##'
##' @param BFmat square Bayes factor matrix or a named vector of posterior probabilities
##' @param step.size how much Bayes factor support is required to accept a more complex model, see Details
##' @param expectedNumberOfShifts expected number of shifts under the prior (only needed for \code{inputType = 'postProb'})
##' @param inputType describes the input: \code{'matrix'} or \code{'postProb'}
##'
##' @details 
##' stepBF takes either a square Bayes factor matrix (such as output by \code{\link{computeBayesFactors}}) or a named 
##' vector of posterior probabilities. If posterior probabilities are supplied, the model prior
##' (\code{expectedNumberOfShifts}) must also be provided.
##' If the input is a Bayes factor matrix, specify \code{inputType = 'matrix'}, otherwise if the input is
##' a named vector of posterior probabilities, specify \code{inputType = 'postProb'}.
##' 
##' The \code{step.size} argument is how much Bayes factor support is needed to accept a more complex model. 
##' By default, this value is 1, so any more complex model that has a better Bayes factor than the previous model 
##' will be accepted. Increasing the step size greatly reduces the Type I error at the cost of inflating Type II 
##' error. So, with increasing step.size, you will infer fewer shifts.
##'
##' @return
##' a list of 3 items: the number of shifts for the best model, the number of shifts for the second best model,
##' and the Bayes factor support for the best model over the second best.
##' @author Jonathan Mitchell
##'
##' @references \url{http://bamm-project.org/}
##' @seealso \code{\link{computeBayesFactors}}
##'
##' @examples
##' data(mcmc.whales)
##' # remove 10% burnin
##' mcmc.whales <- mcmc.whales[floor(0.1 * nrow(mcmc.whales)):nrow(mcmc.whales), ]
##' # from a square matrix of Bayes factor values (inputType = 'matrix')
##' bfmat <- computeBayesFactors(mcmc.whales, expectedNumberOfShifts = 1, burnin = 0)
##' stepBF(bfmat, step.size = 1, inputType = 'matrix')
##' # or from a vector of posterior probabilities (inputType = 'postProb')
##' postProb <- table(mcmc.whales$N_shifts) / nrow(mcmc.whales)
##' stepBF(postProb, step.size = 1, inputType = 'postProb')
##' 
##' @export


stepBF <- function(BFmat, step.size = 20, expectedNumberOfShifts = 1, inputType = 'matrix') {
	inputType <-  match.arg(inputType, c('matrix', 'postProb'))
	if (inherits(BFmat, "table")) {
		BFmat <- setNames(as.vector(BFmat), names(BFmat))
	}
	if (inputType == 'postProb') {
		if (inherits(BFmat, "matrix")) {
			stop("If inputType is 'postProb', please provide a vector of posterior probabilities.")
		}
		post <- BFmat
		prior <- dgeom(as.numeric(names(post)), prob = 1 / (1 + expectedNumberOfShifts))
		BFmat <- matrix(0, nrow = length(post), ncol = length(post))
		rownames(BFmat) <- names(post)
		colnames(BFmat) <- names(post)
		for (i in 1:length(post)) {
			for (j in 1:length(post)) {
				if (post[j] == 0) {
					BFmat[i,j] <- NA
				} else {
					prior_odds <- prior[i] / prior[j]
					post_odds <- post[i] / post[j]
					BFmat[i,j] <- post_odds * (1/prior_odds)
				}		
			}
		}
	}
	if (inherits(BFmat, 'numeric') | !identical(ncol(BFmat), nrow(BFmat))) {
		stop("Bayes factor matrix, BFmat, must be square, with each cell (i,j) representing the BF of model i relative to model j.")
	}
	if (step.size < 1) {
		stop("Step size is less than 1! This means you will move to more complex models even when they aren't supported by any evidence!")
	}
	if (ncol(BFmat) == 1) {
		bestModel <- as.numeric(colnames(BFmat)[1])
		secondModel <- NA
		BF <- NA
	} else if (ncol(BFmat) > 1) {
		Stop <- ncol(BFmat) - 1
		Out <- vector('numeric', length = Stop)
		for (i in 1:Stop) {
			Out[i] <- ifelse(BFmat[i + 1, i] > step.size, BFmat[i + 1, i], 0)
		}
		best <- which(Out == 0)[1]
		second <- best - 1
		bestModel <- as.integer(colnames(BFmat)[best])
		secondModel <- as.integer(colnames(BFmat)[second])
		BF <- as.numeric(BFmat[best, second])
	}
	return(list(bestModel = bestModel, secondModel = secondModel, BF = BF))
}