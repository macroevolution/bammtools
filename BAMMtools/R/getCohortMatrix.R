#############################################################
#
#	getCohortMatrix(....)
#
# 	Each entry of this matrix represents the expected 
#	probability that a pair[i, j] of tips will have the same 
#	rate parameters due to BAMM model.
#
# 	Should modify this to allow exponential, spherical,
#		and other possible correlation structures.
#	Need to make a corStruct class that works with this
#		for GLS analyses

# getCohortMatrix <- function(ephy) {

	# if (!'bammdata' %in% class(ephy)) {
		# stop("Object ephy must be of class bammdata\n");
	# }
	
	# TOL <- 0.0001;
	# corMat <- matrix(0, nrow=length(ephy$tip.label), ncol=length(ephy$tip.label));
	# n <- length(ephy$numberEvents);
	# for (i in 1:length(ephy$tipStates)) {
		# dd <- dist(ephy$tipStates[[i]]);
		# cmat <- as.matrix(dd);	
		# corMat <- corMat + (cmat < TOL)/n;
		
	# }
	# #rownames(corMat) <- ephy$tip.label;
	# #colnames(corMat) <- ephy$tip.label;	
	# dimnames(corMat)[1:2] <- list(ephy$tip.label);
	# #return(corMat/length(ephy$numberEvents));
	# return(corMat);
# }

##' @title Compute the pairwise correlation in rate regimes between all tips
##'     in a \code{bammdata} object
##'
##' @description Takes a \code{bammdata} object and computes the pairwise
##'     correlation in evolutionary rate regimes between all tips in the
##'     phylogeny. This can be used to identify cohorts of taxa that share
##'     common macroevolutionary rate parameters. It can also be used to
##'     construct a correlation matrix for GLS analyses using
##'     \code{BAMM}-estimated tip rates of speciation, extinction, or
##'     phenotypic evolution.
##'
##' @param ephy An object of class \code{bammdata}.
##'
##' @details The cohort matrix is important for interpreting and visualizing
##'     macroevolutionary dynamics. Each entry [i, j] of the cohort matrix is
##'     the probability that taxon i and taxon j share a common
##'     macroevolutionary rate regime. To compute this, we simply tabulate the
##'     percentage of samples from the posterior where taxon i and taxon j
##'     were placed in the same rate regime. If there is no rate heterogeneity
##'     in the dataset (e.g., the data are best explained by a single rate
##'     regime), then all species will tend to share the same rate regime and
##'     all values of the cohort matrix will approach 1. 
##'
##'     A value of 0 between any two taxa means that at least one rate shift
##'     occurred on the nodal path connecting them in 100\% of samples from
##'     the posterior. A value of 0.50 would imply that 50\% of samples from
##'     the posterior included a rate shift on the path connecting taxa i and
##'     j. See below (Examples) for an illustration of this.  
##'
##' @return A numeric matrix of dimension k x k, where k is the number of
##'     species in the phylogeny included in the \code{bammdata} object.
##'     Species names are included as row names and column names. The matrix
##'     is symmetric, such that the values for entry [i , j] will mirror those
##'     for [j , i].
##'
##' @author Dan Rabosky
##'
##' @references \url{http://bamm-project.org/}
##'
##' @examples
##' data(whales, events.whales)
##' ed <- getEventData(whales, events.whales, nsamples=500)
##' 
##' cormat <- getCohortMatrix(ed)
##' 
##' dim(cormat)
##' hist(cormat, breaks=50)
##' @keywords models
##' @export
getCohortMatrix <- function(ephy) {
	if (!inherits(ephy, 'bammdata')) {
		stop("Object ephy must be of class bammdata\n");
	}
	tipStates <- unlist(ephy$tipStates);
	Ntips <- length(ephy$tip.label);
	Nsamples <- length(ephy$tipStates);
	mat <- .C("cohort_matrix",as.integer(tipStates), as.integer(Nsamples), as.integer(Ntips), double(Ntips*Ntips), PACKAGE = "BAMMtools")[[4]];
	dim(mat) <- c(Ntips, Ntips);
	dimnames(mat) <- list(ephy$tip.label, ephy$tip.label);
	diag(mat) <- rep(1.0,Ntips);
	return(mat);
}

