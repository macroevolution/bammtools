##' @title Compute Bayes Factors
##'
##' @description Computes pairwise Bayes factors for a set of
##'     macroevolutionary models sampled using \code{BAMM}, using MCMC
##'     simulation output.
##'
##' @param postdata Filename for the MCMC output file from a \code{BAMM} run.
##'     Alternatively, a dataframe containing this information.
##' @param expectedNumberOfShifts Expected number of shifts under the prior.
##' @param burnin What fraction of samples to discard from postdata as burnin?
##' @param \dots Additional arguments to computeBayesFactors.
##'
##' @details This function returns a matrix of pairwise Bayes factors, where
##'     the Bayes factor is the ratio of marginal likelihoods between two
##'     models M_{i} and M_{j}. Numerator models are given as rows, and
##'     denominator models as columns. Row names and column names give the
##'     number of shifts in the corresponding model. Suppose you have an
##'     output matrix with row and column names 0:3 (0, 1, 2, 3). Model 0 is a
##'     model with just a single process (starting at the root), and no
##'     among-lineage rate heterogeneity. 
##'
##'     If \code{computeBayesFactors} gives a matrix \code{mm}, and
##'     \code{mm[2,1]} is 10.0, this implies Bayes factor evidence of 10 in
##'     favor of the 2nd row model (a model with 1 process; e.g.,
##'     \code{rownames(mm)[2]}) over the first column model (a model with a
##'     single process).
##'
##'     This function will only compute Bayes factors between models which
##'     were actually sampled during simulation of the posterior. Hence, if
##'     a model has such low probability that it is never visited by
##'     \code{BAMM} during the simulation of the posterior, it will be
##'     impossible to estimate its posterior probability (and thus, you will
##'     get no Bayes factors involving this particular model). This is likely
##'     to change in the future with more robust methods for estimating
##'     posterior probabilities in the tails of the distribution. 
##'
##' @return A matrix of pairwise Bayes factors between models.
##'
##' @author Dan Rabosky
##'
##' @examples
##' data(mcmc.whales)
##' computeBayesFactors(mcmc.whales, expectedNumberOfShifts = 1, burnin = 0.1)
##' @keywords models
##' @export
computeBayesFactors <- function(postdata, expectedNumberOfShifts, burnin = 0.1, ...){

	if (hasArg("strict") | hasArg("threshpost") | hasArg("threshprior") | hasArg("nbprior") | hasArg("priordata") | hasArg("modelset")){
 		cat("Error - you have specified some argument names that have been deprecated\n");
 		cat("in this version of BAMMtools. Check the help file on this function\n");
 		cat("to see what has changed\n\n");
		stop();
	}

	if (inherits(postdata, 'character')){
		dpost <- read.csv(postdata, header=T);
	} else if (inherits(postdata, 'data.frame')){
		dpost <- postdata;
	} else{
		stop("invalid postdata argument (wrong class) in computeBayesFactors\n");
	}
 
	dpost <- dpost[floor(burnin*nrow(dpost)):nrow(dpost), ];
 
	tx <- table(dpost$N_shifts) / nrow(dpost);
	
	post <- data.frame(N_shifts=as.numeric(names(tx)), prob=as.numeric(tx));
	
	ux <- as.numeric(names(tx))
	if (length(ux) <= 1){
		cat("Not enough models sampled in simulation of posterior\n")
		cat("You must have valid posterior probabilities for at least 2 models\n")
		cat("to use this function'\n")
	}
	
	pp <- (1 / (1 + expectedNumberOfShifts))
	
	prior <- dgeom(ux, prob = pp)
	names(prior) <- ux
	
 	mm <- matrix(NA, nrow=length(prior), ncol=length(prior));
	rownames(mm) <- names(prior);
	colnames(mm) <- names(prior);
 
	for (i in 1:length(prior)){
			mi <- ux[i];
		for (j in 1:length(prior)){
			
			mj <- ux[j];
				
			prior_odds <- prior[i] / prior[j]
 
			post_odds <- post$prob[post$N_shifts == mi] / post$prob[post$N_shifts == mj];
			
			mm[i,j] <- post_odds * (1 / prior_odds);
		}
	}
	return(mm);
}



