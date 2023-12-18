#############################################################
#
#	getMarginalBranchRateMatrix(....)
#
#	get matrix of marginal rates on each branch for each sample from posterior
# 	
#	This function can handle either a 'bammdata' object or a multiphylo object (i.e., list of trees)

##' @title Compute mean branch rates for \code{bammdata} object
##'
##' @description For each sample in the posterior, computes the mean rate for
##'     each branch in the focal phylogeny (speciation, extinction, trait
##'     evolution). If the \code{bammdata} object contains \emph{nsamples}
##'     samples and the target phylogeny has \emph{nbranches} branches, the
##'     function will compute a matrix of \emph{nbranches} x \emph{nsamples}.
##'
##' @param ephy An object of class \code{bammdata}.
##' @param verbose Print progress during processing of \code{bammdata} object.
##'
##' @details If a \code{type = 'diversification'} \code{bammdata} object is
##'     passed as an argument, the function will return matrices for both
##'     speciation and extinction. If \code{type = 'trait'} object, the matrix
##'     will simply be the corresponding phenotypic rates. Branch-specific
##'     rates are the mean rates computed by integrating the relevant
##'     rate-through-time function along each branch, then dividing by the
##'     length of the branch.
##'
##' @return Returns a list with the following components:
##'     \itemize{
##'         \item lambda_branch_matrix: A \code{nbranches x nsamples} matrix
##'             of mean speciation rates for each branch.
##'         \item mu_branch_matrix: A \code{nbranches x nsamples} matrix of
##'             mean extinction rates for each branch.
##'         \item beta_branch_matrix: A \code{nbranches x nsamples} matrix of
##'             mean phenotypic rates for each branch.
##' }
##'
##' @author Dan Rabosky
##'
##' @examples
##' data(whales)
##' data(events.whales)
##' ed <- getEventData(whales, events.whales, nsamples = 10)
##' mbr <- getMarginalBranchRateMatrix(ed)
##' dim(mbr$lambda_branch_matrix)
##' @keywords models
##' @export
getMarginalBranchRateMatrix <- function(ephy, verbose = FALSE) {
	
	if (!inherits(ephy, 'bammdata')) {
		stop("Object must be of class bammdata\n");
	}

	lammat <- matrix(0, ncol=length(ephy$eventBranchSegs), nrow=nrow(ephy$edge));
	mumat <- matrix(0, ncol=length(ephy$eventBranchSegs), nrow=nrow(ephy$edge));
	
	for (i in 1:length(ephy$eventBranchSegs)) {
		if (verbose) {
			cat('Processing sample ', i, '\n');
		}
		esegs <- ephy$eventBranchSegs[[i]];
		events <- ephy$eventData[[i]];
		events <- events[order(events$index), ];			
		
		# relative start time of each seg, to event:
		relsegmentstart <- esegs[,2] - ephy$eventData[[i]]$time[esegs[,4]];
		relsegmentend <- esegs[,3] - ephy$eventData[[i]]$time[esegs[,4]];
		lam1 <- ephy$eventData[[i]]$lam1[esegs[,4]];
		lam2 <-  ephy$eventData[[i]]$lam2[esegs[,4]];
		mu1 <-  ephy$eventData[[i]]$mu1[esegs[,4]];
		mu2 <-  ephy$eventData[[i]]$mu2[esegs[,4]];
 
		lamint <- timeIntegratedBranchRate(relsegmentstart, relsegmentend, lam1, lam2);
		muint <- timeIntegratedBranchRate(relsegmentstart, relsegmentend, mu1, mu2);
		seglengths <- esegs[,3] - esegs[,2];	
				
		for (k in 1:nrow(ephy$edge)) {
			isRightBranch <- esegs[,1] == ephy$edge[k,2];
			lammat[k, i] <- sum(lamint[isRightBranch]) / sum(seglengths[isRightBranch]);
			mumat[k, i] <- sum(muint[isRightBranch]) / sum(seglengths[isRightBranch]);
		}
	}
	
	if (ephy$type == 'diversification') {
		return(list(lambda_branch_matrix = lammat, mu_branch_matrix = mumat));
	}
	if (ephy$type == 'trait') {
		return(list(beta_branch_matrix = lammat));
	}
	  
}
