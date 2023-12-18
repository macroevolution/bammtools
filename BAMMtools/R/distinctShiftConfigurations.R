# Drop all nodes from event data with marginal probs < 0.05
# test is unique
#  if so, add to list
# Returns:
#	$marg.probs = marginal probs for nodes
#	$marginal_odds_ratio = branch-specific (marginal) posterior:prior odds ratios associated with 1 or more shifts
#	$shifts = unique shift sets
#	$samplesets = list of sample indices that reduce to each of the unique shift sets
#	$frequency = vector of frequencies of each shift configuration
#	$threshold =  (marginal) posterior:prior odds ratio threshold for shifts
#	
#	Results are sorted by frequency. 
#	$frequency[1] gives the most common shift config sampled
#	$shifts[[1]] gives the corresponding node indices for that configuration
#	$samplesets[[1]] gives the indices of samples with this configuration


##' @title Identify distinct rate shift configurations
##'
##' @description Identify topologically distinct rate shift configurations
##'     that were sampled with \code{BAMM}, and assign each sample in the
##'     posterior to one of the distinct shift configurations.
##'
##' @param ephy An object of class \code{bammdata}.
##' @param expectedNumberOfShifts The expected number of rate shifts under the
##'     prior.
##' @param threshold Threshold value for marginal posterior-to-prior odds
##'     ratios, used to identify branches with elevated shift probabilities
##'     relative to prior (core vs non-core shifts).
##' @param \dots Other arguments to distinctShiftConfigurations (possibly
##'     deprecated args).
##'
##' @details See Rabosky et al (2014) and the \code{BAMM} project website for
##'     details on the nature of distinct shift configurations and especially
##'     the distinction between "core" and "non-core" rate shifts. Note that
##'     branches with elevated marginal posterior probabilities relative to
##'     the prior (marginal odds ratios) cannot be claimed to have
##'     "significant" evidence for a rate shift on the basis of this evidence
##'     alone.  
##'
##' @return An object of class \code{bammshifts}. This is a list with the
##'     following components:
##'     \itemize{
##'         \item marg.probs: A list of the marginal probability of a shift
##'             occurring at each node of the phylogeny for each distinct rate
##'             shift configuration.
##'         \item marginal_odd_ratio: Marginal posterior-to-prior odds ratios
##'             for one or more rate shifts an a given branch.
##'         \item shifts: A list of the set of shift nodes for each distinct
##'             rate configuration.
##'         \item samplesets: A list of sample indices that reduce to each of
##'             the unique shift sets.
##'         \item frequency: A vector of frequencies of each distinct shift
##'             configuration.
##'         \item coreshifts: A vector of node numbers corresponding to the
##'             core shifts. All of these nodes have a marginal odds ratio of
##'             at least \code{threshold} supporting a rate shift.
##'         \item threshold: A single numeric value giving the marginal
##'             posterior:prior odds ratio threshold used during enumeration
##'             of distinct shift configurations.
##'     }
##'     Results are sorted by frequency:
##'
##'     $frequency[1] gives the most common shift configuration sampled.
##'
##'     $shifts[[1]] gives the corresponding node indices for that
##'     configuration.
##'
##'     $samplesets[[1]] gives the indices of samples with this configuration.
##'
##' @author Dan Rabosky
##'
##' @seealso \code{\link{plot.bammshifts}}, \code{\link{credibleShiftSet}}
##'
##' @examples
##' data(whales, events.whales)
##' ed <- getEventData(whales, events.whales, burnin=0.25, nsamples=500)
##' 
##' sc <- distinctShiftConfigurations(ed, expectedNumberOfShifts = 1,
##'                                   threshold = 5)
##' 
##' plot(sc, ed, rank=1)
##' @export
distinctShiftConfigurations <- function(ephy, expectedNumberOfShifts, threshold, ... ) {
  
	or <- marginalOddsRatioBranches(ephy, expectedNumberOfShifts)
	
	mm <- marginalShiftProbsTree(ephy);

	goodnodes <- or$edge[,2][or$edge.length >= threshold];

	xlist <- list();
	for (i in 1:length(ephy$eventData)) {
		xlist[[i]] <- intersect(goodnodes, ephy$eventData[[i]]$node);
	}

	ulist <- list();
	treesets <- list();
	
	ulist[[1]] <- xlist[[1]];
	treesets[[1]] <- 1;
	
	for (i in 2:length(xlist)) {
		lx <- length(ulist);
		#cat(lx, '\n')
		for (k in 1:lx) {
			if (areShiftSetsEqual(ulist[[k]], xlist[[i]])){
				treesets[[k]] <- c(treesets[[k]], i);
				break;	
			} else {
				if (k == length(ulist)){
					xlen <- length(ulist);
					ulist[[xlen + 1]] <- xlist[[i]];
					treesets[[xlen + 1]] <- i;
				}
			}
		}
	}
	
	freqs <- unlist(lapply(treesets, length));
	freqs <- freqs / sum(freqs);
	
	ord <- order(freqs, decreasing=TRUE);
	
	obj <- list();
	obj$marg.probs <- mm$edge.length;  
	names(obj$marg.probs) <- mm$edge[,2]; 
	obj$marginal_odds_ratio <- or$edge.length;
	names(obj$marginal_odds_ratio) <- or$edge[,2];
	obj$shifts <- ulist[ord]; 
	obj$samplesets <- treesets[ord];
	obj$frequency <- freqs[ord];
	obj$coreshifts <- goodnodes;
	obj$threshold <- threshold;

	class(obj) <- 'bammshifts';
	
	return(obj);
}

