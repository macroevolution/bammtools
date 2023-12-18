#############################################################
#
#	getMeanBranchLengthTree(....)
#

##' @title Compute phylogeny with branch lengths equal to corresponding
##'     macroevolutionary rate estimates
##'
##' @description Takes a \code{bammdata} object and computes a phylogenetic
##'     tree where branch lengths are equal to the mean of the marginal
##'     distributions of rates on each branch. This tree can be plotted to
##'     visualize rate variation across a phylogeny.
##'
##' @param ephy An object of class \code{bammdata}.
##' @param rate The type of rate-tree to be computed. Options: "speciation"
##'     (default), "extinction", "ndr" (net diversification), and "trait".
##' 
##' @return A list with the following components:
##'     \itemize{
##'         \item phy: A phylogenetic tree, topologically identical to the
##'             model tree, but with branch lengths replaced by the mean
##'             (marginal) rates on each branch as estimated from the
##'             posterior samples in the \code{bammdata} object.
##'         \item mean: The mean rate over all branches.
##'         \item median: the median rate over all branches.
##'     }
##'
##' @author Dan Rabosky
##'
##' @seealso \code{\link{plot.bammdata}}
##'
##' @references \url{http://bamm-project.org/}
##'
##' @examples
##' data(whales)
##' data(events.whales)
##' ed <- getEventData(whales, events.whales, burnin=0.1, nsamples=500)
##' ed2 <- subsetEventData(ed, index = 1:20)
##' ratetree <- getMeanBranchLengthTree(ed2, rate='speciation')
##' plot(ratetree$phy, show.tip.label=FALSE)
##' @keywords models
##' @export
getMeanBranchLengthTree <- function(ephy, rate='speciation') {
	
	if (inherits(ephy, 'bammdata')) {
		v <- as.phylo.bammdata(ephy);
	}
	obj <- getMarginalBranchRateMatrix(ephy,verbose=FALSE);

	if (ephy$type == 'diversification'){
		if (rate == 'speciation'){
			el <- rowMeans(obj$lambda_branch_matrix);
		}else if (rate == 'extinction'){
			el <- rowMeans(obj$mu_branch_matrix);
		}else if (rate == 'ndr'){
			el  <- rowMeans(obj$lambda_branch_matrix) - rowMeans(obj$mu_branch_matrix);				
		}else{
			stop("invalid rate specification in getMeanBranchLengthTree");
		}
		
	}else if (ephy$type == 'trait'){
		el <- rowMeans(obj$beta_branch_matrix);
	
	}else{
		stop("error in getMeanBranchLengthTree - \nproblem with supplied ephy object");
	}
	
	v$edge.length <- el;
	tmp <- list();
	tmp$phy <- v;
	tmp$median <- median(el);
	tmp$mean <- mean(el);		
		
	
	return(tmp);		
}
