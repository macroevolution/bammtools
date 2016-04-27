
# returns a phylogenetic tree where 
# branch lengths are equal to the marginal 
# odds ratio (posterior : prior) for a given branch. 
# It is marginal in the sense that it is not independent 
# of values for other branches.

##' @title Ratio of (marginal) posterior-to-prior probabilities on individual
##'     branches
##'
##' @description Compute marginal posterior-to-prior odds ratio associated
##'     with observing one or more rate shift on a given branch.
##'
##' @param ephy An object of class \code{bammdata}.
##' @param expectedNumberOfShifts Expected number of shifts under the prior
##'     alone.
##'
##' @details This function returns a copy of a phylogenetic tree where each
##'     branch length is equal to the marginal odds ratio in favor of a rate
##'     shift on a particular branch. These cannot be interpreted as evidence
##'     for a rate shift in an absolute sense. As explained on the website,
##'     they are a marginal odds ratio. This function is provided primarily
##'     for the purpose of distinguishing core and non-core shifts.
##'
##' @return A object of class \code{phylo} but where each branch length is
##'     equal to the marginal shift odds on each branch.
##'
##' @author Dan Rabosky
##'
##' @seealso \code{\link{getBranchShiftPriors}},
##'     \code{\link{distinctShiftConfigurations}},
##'     \code{\link{credibleShiftSet}}
##'
##' @examples
##' #Produce a blank template control file
##' generateControlFile(file = 'traitcontrol.txt', type='trait')
##'
##' #Produce a customized control file
##' data(whales)
##'
##' #get bamm priors to supply to control file
##' priors <- setBAMMpriors(whales, outfile = NULL)
##'
##' generateControlFile(file = 'divcontrol.txt', params = list(
##' 	treefile = 'whales.tre',
##' 	globalSamplingFraction = '1',
##' 	numberOfGenerations = '100000',
##' 	overwrite = '1',
##' 	lambdaInitPrior = as.numeric(priors['lambdaInitPrior']),
##' 	lambdaShiftPrior = as.numeric(priors['lambdaShiftPrior']),
##' 	muInitPrior = as.numeric(priors['muInitPrior']),
##' 	expectedNumberOfShifts = '1'))
##' @export
marginalOddsRatioBranches <- function(ephy, expectedNumberOfShifts) {
	
	tree_post <- marginalShiftProbsTree(ephy);
	tree_prior <- getBranchShiftPriors(as.phylo.bammdata(ephy), expectedNumberOfShifts);
	
	post_shift <- tree_post$edge.length;
	prior_shift <- tree_prior$edge.length;
 
	oddsratio <- post_shift / prior_shift;

	tree_post$edge.length <- oddsratio;
	return(tree_post);	
	
}



