#############################################################
#
#	getShiftNodesFromIndex (....)
#
#	Args: ephy	=	object of class 'bammdata'
#		  index =   the index of the sample you wish to view, e.g.,
#					 if index = 5, this will give you the nodes subtending
#                    all branches with rate shifts for the 5th sample
#					 from the posterior in your 'bammdata' object.								 
#
#	
#	Returns: 		- a vector of the nodes where shifts occurred, excluding the root.
#					Note: if NO shifts occured, this will return a 
#							numeric object of length zero
# 

##' @title Identify nodes associated with rate shifts from \code{bammdata}
##'     object
##'
##' @description Find the node numbers associated with rate shifts for a
##'     specified sample from the posterior distribution contained in a
##'     \code{bammdata} object.
##'
##' @param ephy A \code{bammdata} object.
##' @param index The index value of the posterior sample from which you want
##'     to identify shift nodes. This is \emph{not} the same as the actual
##'     generation number of the MCMC sample. If your \code{bammdata} object
##'     contains 100 samples from the posterior distribution, the value of
##'     \code{index} must range from 1 to 100.
##'
##' @return A vector of nodes (excluding the root) that define branches on
##'     which shifts occurred for the specified sample from the posterior.
##'     Will return a numeric of length 0 if no non-root shifts occur in the
##'     specified sample.
##'
##' @author Dan Rabosky
##'
##' @seealso \code{\link{addBAMMshifts}}, \code{\link{plot.bammdata}},
##' \code{\link{maximumShiftCredibility}}
##'
##' @examples
##' data(whales, events.whales)
##' ed <- getEventData(whales, events.whales, burnin=0.1, nsamples=500)
##' 
##' # Get the maximum shift credibility configuration:
##' msc <- maximumShiftCredibility(ed)
##' 
##' # Get the nodes at which shifts occurred in the 
##' # maximum shift credibility configuration:
##' 
##' getShiftNodesFromIndex(ed, index=msc$sampleindex)
##' @keywords models
##' @export
getShiftNodesFromIndex <- function(ephy, index) {
	
	if (index > length(ephy$eventData)) {
		cat("Error. Attempting to access non-existent element from 'bammdata' object\n");
		cat("You have << ", length(ephy$eventData), " >>> samples in your 'bammdata' object\n");
		cat("Value of index must be no greater than the number of samples\n\n");
		stop();
	}
	root <- length(ephy$tip.label) + 1;
	nodes <- ephy$eventData[[index]]$node;
	nodes <- nodes[nodes != root];

	return(nodes);
}
