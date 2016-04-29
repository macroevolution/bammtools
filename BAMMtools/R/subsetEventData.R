### subsetEventData

##' @title Subset a \code{bammdata} object
##'
##' @description Subsets a \code{bammdata} object. Returns a \code{bammdata}
##'     object after extracting a specified set of samples from the posterior.
##'
##' @param ephy An object of class \code{bammdata}.
##' @param index A vector of integers corresponding to samples to be extracted
##'     from the posterior distribution of shift configurations included in
##'     the \code{bammdata} object.
##'
##' @details This will result in an error if you attempt to access samples
##'     that do not exist in the \code{ephy} data object. For example, if your
##'     \code{bammdata} object includes 100 samples from a posterior
##'     distribution sampled with \code{BAMM}, you can only attempt to subset
##'     with index values 1:100.
##'
##' @author Dan Rabosky
##'
##' @seealso \code{\link{plot.bammdata}}, \code{\link{getCohortMatrix}},
##'     \code{\link{image}}
##'
##' @examples
##' data(whales, events.whales)
##' ed <- getEventData(whales, events.whales, nsamples=500)
##' ed2 <- subsetEventData(ed, index=1)
##' plot(ed2) 
##' addBAMMshifts(ed2, cex=2)
##' @keywords manip
##' @export
subsetEventData <- function(ephy, index) {
	
	if (class(ephy) != 'bammdata') {
		stop("Object ephy must be of class bammdata\n");
	}
	
	nsamples <- length(ephy$eventData);

	ss <- which((1:nsamples) %in% index);
	badsubset <- setdiff(index, 1:nsamples);
	
	if (length(badsubset) > 0) {
		cat("Bad call to BAMMtools::subsetEventData\n");
		cat("You have << ", nsamples, " >> samples in your data object \n");
		stop("Attempt to access invalid samples. Check index.")
	}
	
	obj <- list();
	obj$edge <- ephy$edge;
	obj$Nnode <- ephy$Nnode;
	obj$tip.label <- ephy$tip.label;
	obj$edge.length <- ephy$edge.length;
	obj$begin <- ephy$begin;
	obj$end <- ephy$end;
	obj$downseq <- ephy$downseq;
	obj$lastvisit <- ephy$lastvisit;
	
	obj$numberEvents <- ephy$numberEvents[ss];
	obj$eventData <- ephy$eventData[ss];	
	obj$eventVectors <- ephy$eventVectors[ss];
	obj$tipStates <- ephy$tipStates[ss]
	obj$tipLambda <- ephy$tipLambda[ss];
	obj$tipMu <- ephy$tipMu[ss];
	obj$eventBranchSegs <- ephy$eventBranchSegs[ss];
	
	obj$meanTipLambda <- ephy$meanTipLambda;
	obj$meanTipMu <- ephy$meanTipMu;
	
	obj$type <- ephy$type;
	
	class(obj) <- 'bammdata';
	attributes(obj)$order = attributes(ephy)$order;	
	return(obj);
}
