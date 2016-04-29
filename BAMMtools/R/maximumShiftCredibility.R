#############################################################
#
#	maximumShiftCredibility(....)
#
#	Args: ephy	=	object of class 'bammdata'
#	
#
#	maximize		=	'sum', = sum of branch probabilities for each tree
#						'product', = product of branch probabilities for each tree
#	
#	Returns: 		- bestconfigs: a list of length equal the number of 
#						unique shift configurations in the maximum shift
#						credibility set. Each element is a vector of sample
#						indices from the 'bammdata' object with identical
#						shift configurations.
#					  
#					- A vector of optimality scores for all other samples 
#						in posterior from the 'bammdata' object.
#
#					- sampleindex: a representative index for samples from each 
#						set of unique shift configurations. The length of this vector
#						is equal to the length of the bestconfigs list. If this vector was
#						sampleindex = c(2, 20, 50), this would mean that there are 3 distinct
#						sets of shift configurations with equal credibility under the optimality 
#						criterion. More commonly, a single shift configuration will be dominant, and 
#						although the length of bestconfigs[[1]] may be greater than 1, the sampleindex
#						vector will contain a single representative event from that set.
#
#	See example file.
#   This is analogous to the maximum clade credibility tree from a 
#		Bayesian phylogenetic analysis.

##' @title Estimate maximum shift credibility configuration
##'
##' @description This is one estimate of the "best" rate shift configuration,
##'     considering only those shift configurations that were actually sampled
##'     using \code{BAMM}'s reversible jump MCMC simulator. This is analogous
##'     to the "maximum clade credibility tree" from a Bayesian phylogenetic
##'     analysis. It is not necessarily the same as the shift configuration
##'     with the maximum a posteriori probability.
##'
##' @param ephy An object of class \code{BAMMdata}.
##' @param maximize Maximize the marginal probability of the product or sum of
##'     branch-specific shifts.
##'
##' @details This is one point estimate of the overall "best" rate shift
##'     configuration. Following an MCMC simulation, the marginal shift
##'     probabilities on each individual branch are computed using
##'     \code{\link{marginalShiftProbsTree}}. The shift configuration that
##'     maximizes the product (or sum, if specified) of these marginal
##'     branch-specific shift probabilities is the \emph{maximum shift
##'     credibility configuration}.  
##' 
##'     This option is only recommended if you have no clear "winner" in your
##'     credible set of shift configurations (see
##'     \code{\link{credibleShiftSet}}). If you have a number of
##'     largely-equiprobable shift configurations in your 95\% credible set,
##'     you may wish to try this function as an alternative for identifying a
##'     single best shift configuration. Otherwise, it is recommended that you
##'     present the shift configuration with the maximum a posteriori
##'     probability (see \code{\link{getBestShiftConfiguration}}).
##'
##' @return A list with the following components:
##'     \itemize{
##'         \item{bestconfigs} {A vector of the index values of MCMC samples
##'             with shift configurations equal to the maximum. Usually, more
##'             than one state sampled during the MCMC simulation will have an
##'             identical (maximized) marginal probability. All samples given
##'             in this vector will have an identical shift configuration.}
##'         \item{scores} {The optimality score (product or sum of marginal
##'             shift probabilities) for all sampled shift configurations in
##'             the \code{BAMMdata} object.}
##'         \item{optimalityType} {Whether the product or sum of marginal
##'             shift probabilities was used to compute the maximum shift
##'             credibility configuration.}
##'         \item{sampleindex} {A representative sample that is equal to the
##'             maximum shift credibility configuration (e.g., this can be
##'             plotted with \code{\link{addBAMMshifts}}).}
##'     }
##'
##' @author Dan Rabosky
##'
##' @seealso \code{\link{marginalShiftProbsTree}},
##'     \code{\link{addBAMMshifts}}, \code{\link{cumulativeShiftProbsTree}},
##'     \code{\link{credibleShiftSet}},
##'     \code{\link{getBestShiftConfiguration}}
##'
##' @examples
##' data(whales, events.whales)
##' ed <- getEventData(whales, events.whales, burnin=0.25, nsamples=500)
##' best_config <- maximumShiftCredibility(ed)
##' plot(ed)
##' addBAMMshifts(ed, method='phylogram', index=best_config$sampleindex)
##' @keywords manip graphics
##' @export
maximumShiftCredibility <- function(ephy, maximize = 'product') {

	if (!'bammdata' %in% class(ephy)) {
		stop("Object ephy must be of class bammdata\n");
	}			
	
	probvec <- numeric(length(ephy$eventData));
	
	mtree <- marginalShiftProbsTree(ephy);
	
	#mtree$edge.length[mtree$edge.length < threshold] <- 0;
	
	px <- mtree$edge.length;
	
	ttx <- table(ephy$numberEvents) / length(ephy$numberEvents);
	
	for (i in 1:length(ephy$eventData)) {
		
		# posterior probabilities here:
		proc_prob <- ttx[as.character(ephy$numberEvents[i])]; 
		
		hasShift <- ephy$edge[,2] %in% ephy$eventData[[i]]$node;
		branchprobs <- (hasShift)*px  + (!hasShift)*(1 - px) ;
		if (maximize == 'product') {
			probvec[i] <- log(proc_prob) + sum(log(branchprobs));
		} else if (maximize == 'sum') {
			probvec[i] <- proc_prob * sum(branchprobs);
		} else {
			stop("Unsupported optimize criterion in maximumShiftCredibilityTree");
		}
	}
	
	best <- which(probvec == max(probvec));
	
	# Now test for multiple trees with same log-prob:
	bestconfigs <- list();
		
	index <- 0;	
	while (length(best) > 0) {
		index <- index + 1;	
		lv <- logical(length = length(best));
		for (i in 1:length(best)) {
			lv[i] <- areEventConfigurationsIdentical(ephy, best[1], best[i]);
		}
		bestconfigs[[index]] <- best[lv];
		best <- best[!lv];
	}
	
	sampleindex <- numeric(length(bestconfigs));
	for (i in 1:length(bestconfigs)) {
		sampleindex[i] <- bestconfigs[[i]][1];
	}
	
	obj <- list();
	obj$bestconfigs <- bestconfigs;
	obj$scores <- probvec;
	obj$optimalityType = maximize;
 	obj$sampleindex <- sampleindex;
	return(obj);
}
