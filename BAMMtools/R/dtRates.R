############################################
#	dtRates(ephy,tau)
#
#	A function to calculate approximations of
#	mean instantaneous speciation/extinction rates or
#	phenotypic rates along each branch.
#
#	Arguments: ephy = a bammdata object
#	           tau = fraction of tree height for approximation (e.g. 0.01).
#	                 This is the step size over which rates are calculated along
#	                 a branch, 0.01 corresponds to a step size of 1% of tree height.
#	           ism = index of posterior sample(s). Currently may be NULL or 
#	                 a vector of integer values.  if NULL the function will use all 
#	                 posterior samples, otherwise it will use only
#	                 the samples corresponding to the indices in ism,
#	                 e.g. 50, e.g. 50:100.
#
#	Returns: an ephy object with a list appended containing a vector of branch
#			 rates and the step size used for calculation.

##' @title Calculate macroevolutionary rate changes on a phylogeny from
##'     \code{BAMM} output
##'
##' @description \code{dtRates} calculates the mean of the marginal posterior
##'     density of the rate of speciation/extinction or trait evolution for
##'     small segments along each branch in a phylogeny.
##'
##' @param ephy An object of class \code{bammdata}.
##' @param tau A numeric that specifies the size (as a fraction of tree
##'     height) of the segments that each branch will be discretized into.
##' @param ism An integer vector indexing which posterior samples to include
##'     in the calculation.
##' @param tmat A logical. If \code{TRUE} the matrix of branch segments is
##'     returned.
##'
##' @details \code{dtRates} bins the phylogeny into windows of time and
##'     calculates average rates of speciation/extinction or phenotypic
##'     evolution along each segment of a branch within a window. The width of
##'     each window is determined by \code{tau}. \code{tau} is a fraction of
##'     the root to tip distance so a value of \code{tau = 0.01} bins the
##'     phylogeny into 100 time windows of equal width. 
##'
##' @return A \code{bammdata} object with a new component named "dtrates", which
##'     is a list with two or three components:
##'     \itemize{
##'         \item{tau} {The parameter value of \code{tau} used in the
##'             calculation.}
##'         \item{rates} {If \code{ephy$type = "trait"}: a numeric vector with
##'             the phenotypic rates of each segment on each branch. If
##'             \code{ephy$type = "diversification"}: a list with two
##'             components. The first component is a numeric vector of
##'             speciation rates. The second component is a numeric vector of
##'             extinction rates.}
##'         \item{tmat} {A matrix of the starting and ending times of the
##'             segments on each branch. Only if \code{tmat = TRUE}.}
##' }
##'
##' @note If there are zero length branches in the input tree \code{NA}s will
##'     result.
##'
##' @author Mike Grundler
##'
##' @seealso \code{\link{plot.bammdata}}
##'
##' @references \url{http://bamm-project.org/}
##'
##' @examples
##' data(whales, events.whales)
##' ed <- getEventData(whales, events.whales, burnin=0.25, nsamples=500)
##' 
##' # use all posterior samples
##' ed <- dtRates(ed, tau=0.01)
##' 
##' # use specified range of posterior samples
##' ed <- dtRates(ed, tau=0.01, ism=50:150)
##' @keywords graphics
##' @export
dtRates <- function (ephy, tau, ism = NULL, tmat = FALSE) {
    if (!inherits(ephy, "bammdata")) {
        stop("Object ephy must be of class bammdata");
    }
    if (attributes(ephy)$order != "cladewise") {
    	stop("Function requires tree in 'cladewise' order");
    }
    ephy$eventBranchSegs <- lapply(ephy$eventBranchSegs, function(x) x[order(x[,1]), ]);
#   phy <- as.phylo.bammdata(ephy);
#   phy <- getStartStopTimes(phy);   
#   if (is.ultrametric(phy))
#	    tH <- max(branching.times(phy))
#	else
#		tH <- max(NU.branching.times(phy));

	tH <- max(ephy$end);
    segmat <- segMap(ephy, tau);
    #tol = max(1 * 10^-decimals(ephy$eventBranchSegs[[1]][1, 2]),1 * 10^-decimals(ephy$eventBranchSegs[[1]][1, 3]));
    tol <- 0.00001;
    if (storage.mode(segmat) != "double") stop("Error: wrong storage mode in foreign function call");
    if (storage.mode(tol) != "double") stop("Error: wrong storage mode in foreign function call");
    if (storage.mode(ephy) != "list") stop("Error: wrong storage mode in foreign function call");
    if (is.null(ism)) {
        ism <- as.integer(1:length(ephy$eventBranchSegs));
    }
    else ism <- as.integer(ism);
    if (ism[length(ism)] > length(ephy$eventBranchSegs)) {
        warning("Sample index out of range");
        ism <- as.integer(1:length(ephy$eventBranchSegs));
    }
    index <- 1:nrow(segmat);
    rownames(segmat) <- index;
    segmat <- segmat[order(segmat[, 1]), ];
    if (ephy$type == "diversification") {
        dtrates <- .Call("dtrates", ephy, segmat, tol, ism, 0L, PACKAGE = "BAMMtools");
        for (i in 1:2) {
            names(dtrates[[i]]) <- rownames(segmat);
            dtrates[[i]] <- dtrates[[i]][as.character(index)];
            names(dtrates[[i]]) <- NULL;
        }
        if (sum(is.na(dtrates[[1]]))) {
            warning(sprintf("Found %d NA speciation rates. Coercing to zero.", sum(is.na(dtrates[[1]]))));
            dtrates[[1]][is.na(dtrates[[1]])] <- 0;
        }
        if (sum(is.na(dtrates[[2]]))) {
            warning(sprintf("Found %d NA extinction rates. Coercing to zero.", sum(is.na(dtrates[[2]]))));
            dtrates[[2]][is.na(dtrates[[2]])] <- 0;
        }
    }
    else if (ephy$type == "trait") {
        dtrates <- .Call("dtrates", ephy, segmat, tol, ism, 1L, PACKAGE = "BAMMtools");
        names(dtrates) <- rownames(segmat);
        dtrates <- dtrates[as.character(index)];
        names(dtrates) <- NULL;
        if (sum(is.na(dtrates))) {
            warning(sprintf("Found %d NA phenotypic rates. Coercing to zero.", sum(is.na(dtrates))));
            dtrates[is.na(dtrates)] <- 0;
        }
    }
    else {
    	stop("Unrecognized model type");
    }
	if (tmat) {
		segmat <- segmat[as.character(index),];
		ephy$dtrates <- list(tau = tau, rates = dtrates, tmat = segmat);
		return(ephy);
	}
	ephy$dtrates <- list(tau = tau, rates = dtrates);
	return(ephy);
}


# dtRates = function(ephy, tau, ism = NULL) {
	# if (!'bammdata' %in% class(ephy)) {
		# stop('Object ephy must be of class bammdata');
	# }
	
	# ephy$eventBranchSegs = lapply(ephy$eventBranchSegs, function(x) x[order(x[,1]), ]); 
	
	# phy = as.phylo.bammdata(ephy);
	# phy = getStartStopTimes(phy);
	# #if (attributes(phy)$order != 'cladewise') {
	# #	phy = reorder(phy,'cladewise');
	# #}
	# tH = max(branching.times(phy));
	
	# segmat = segMap(phy$edge[,2],phy$begin/tH,phy$end/tH,tau);
	# segmat[,2] = segmat[,2] * tH;
	# segmat[,3] = segmat[,3] * tH;
	
	# tol = max(1*10^-decimals(ephy$eventBranchSegs[[1]][1,2]),1*10^-decimals(ephy$eventBranchSegs[[1]][1,3]));
	
	# if (storage.mode(segmat) != "double") stop('Exiting');
	# if (storage.mode(tol) != "double") stop('Exiting');
	# if (storage.mode(ephy) != "list") stop('Exiting');
	
	# if (is.null(ism)) ism = as.integer(1:length(ephy$eventBranchSegs)) else ism = as.integer(ism);
	# if (ism[length(ism)] > length(ephy$eventBranchSegs)) {
		# warning("Sample index out of range");
		# ism = as.integer(1:length(ephy$eventBranchSegs));
	# }
	
	# index = 1:nrow(segmat)
	# rownames(segmat) = index;
	# segmat = segmat[order(segmat[,1]),];
	# if (ephy$type == "diversification") {
		# dtrates = .Call("dtrates", ephy, segmat, tol, ism, 0L, PACKAGE="BAMMtools");
		# for (i in 1:2) {
			# names(dtrates[[i]]) = rownames(segmat);
			# dtrates[[i]] = dtrates[[i]][as.character(index)];
			# names(dtrates[[i]]) = NULL;
		# }
	# }
	# else {
		# dtrates = .Call("dtrates", ephy, segmat, tol, ism, 1L, PACKAGE="BAMMtools");
		# names(dtrates) = rownames(segmat);
		# dtrates = dtrates[as.character(index)];
		# names(dtrates) = NULL;
	# }
	# ephy$dtrates = list(tau = tau, rates = dtrates);
	# return(ephy);
# }
