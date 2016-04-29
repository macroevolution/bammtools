##' @title Compute prior odds of a rate shift on each branch of a phylogeny
##'     from BAMM output
##'
##' @description Computes the prior probability of a rate shift event for each
##'     branch. These results are important for identifying topological rate
##'     shift locations on phylogenies with marginal probabilities that exceed
##'     those predicted under the prior alone.
##'
##' @param phy An object of class \code{phylo}.
##' @param expectedNumberOfShifts The expected number of shifts under the
##'     prior.
##'
##' @details This function computes the prior odds on the distribution of
##'     numbers of rate shift events per branch under the prior. It returns an
##'     object which is nothing more than a copy of the original phylogenetic
##'     tree but where each branch length has been replaced by the prior
##'     probability of a rate shift on each branch.
##'
##'     The significance of this function is that it lets us explicitly
##'     determine which branches have shift probabilities that are elevated
##'     relative to the prior expectation.  
##'
##' @return A class \code{phylo} with all the components of the original class
##'     \code{phylo} object, with the following changes:
##'
##'     \item{edge.length}{Branch lengths now represent the prior probability
##'     of a rate shift on each branch.}
##'
##' @author Dan Rabosky
##'
##' @seealso \code{\link{distinctShiftConfigurations}},
##'     \code{\link{plot.bammshifts}}, \code{\link{summary.credibleshiftset}},
##'     \code{\link{plot.credibleshiftset}}, \code{\link{credibleShiftSet}}
##'
##' @references \url{http://bamm-project.org}
##'
##' @examples
##' data(whales)
##' prior_tree1 <- getBranchShiftPriors(whales, expectedNumberOfShifts = 1)
##' prior_tree10 <- getBranchShiftPriors(whales, expectedNumberOfShifts = 10)
##' # plot prior expectations for branches based on these two counts:
##' plot(prior_tree1$edge.length ~ prior_tree10$edge.length, xlim=c(0,0.05),
##'      ylim=c(0,0.05), asp=1)
##' lines(x=c(0,1), y=c(0,1))
##' @keywords models
##' @export
getBranchShiftPriors <- function(phy, expectedNumberOfShifts) {
	
	Nmax <- 1000;
	
	geom_p <- 1 / (expectedNumberOfShifts + 1);
	prior <- dgeom(1:Nmax, geom_p);
 
	pvec <- phy$edge.length / sum(phy$edge.length);
	
	pp <- numeric(length(phy$edge.length));
 
	for (i in 1:length(prior)){
		# probability of getting 0 shifts on branch given ns total 
		#  given the branch lengths etc
		#	weighted by the probability of that shift category
		pp <- pp + (1 - dbinom(0, i, prob=pvec)) * prior[i];
	}
	
	obj <- phy;	
	obj$edge.length <- pp;
	return(obj);
}

