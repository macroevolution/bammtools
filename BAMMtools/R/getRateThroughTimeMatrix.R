#############################################################
#
#	getRateThroughTimeMatrix(....)
#	
#	include or exclude options:
#	start.time		=	 start time (in units before present)
#						 if NULL, starts at root
#	end.time		=	 end time 
#						 if NULL, ends at present
#	nslices			=	 number of time cuts
#	Return
#	list with components:
#					$times	= the time vector
#					$lambda = speciation rate matrix
#					$mu 	= extinction rate matrix
# 					$type   = diversification or trait (needs to be extended to trait data)		
# returns object of class bamm-ratematrix
#	 

##' @title Generate rate-through-time matrix from \code{bammdata} object
##'
##' @description Computes a matrix of macroevolutionary rates at specified
##'     timepoints from a \code{bammdata} object. These rates can be used for
##'     plotting speciation rates (and other rates) through time.
##'
##' @param ephy An object of class \code{bammdata}.
##' @param start.time The start time (in units before present) for the time.
##'     sequence over which rates should be computed. If \code{NULL}, starts
##'     at the root.
##' @param end.time The end time (in units before present) for the time
##'     sequence over which rates should be computed. If \code{NULL}, ends in
##'     the present.
##' @param nslices The number of time points at which to compute rate
##'     estimates (between \code{start.time} and \code{end.time}).
##' @param node Allows user to extract rate-through-time information for the
##'     subtree descended from a specific node. Alternatively, a specified
##'     subtree can be excluded from the rate matrix calculations.
##' @param nodetype Two options: "include" and "exclude". If "include",
##'     computes rate matrix only for the descendants of subtree defined by
##'     \code{node}. If "exclude", computes rate matrix for all background
##'     lineages in tree after excluding the subtree defined by \code{node}.
##'     Ignored if \code{node = NULL}.
##'
##' @details Computes evolutionary rates for each sample in the posterior
##'     included as part of the \code{bammdata} object. Rates are computed by
##'     draping an imaginary grid over the phylogeny, where the grid begins at
##'     \code{start.time} and ends at \code{end.time}, with \code{nslices}
##'     vertical lines through the phylogeny. The mean rate at each point in
##'     time (for a given sample from the posterior) is simply the mean rate
##'     at that time for all branches that are intersected by the grid (see
##'     the grid plot in the examples section).
##'
##'     This function is used by \link{plotRateThroughTime}, but the user can
##'     work directly with the \code{bamm-ratematrix} object for greater
##'     control in plotting rate-through-time trajectories for individual
##'     clades. See \code{examples} for an example of how this can be used to
##'     plot confidence intervals on a rate trajectory using shaded polygons.
##'
##'     The \code{node} options are particularly useful. If you have run
##'     \code{BAMM} on a large phylogeny, you can easily generate the
##'     rate-through-time data for a particular subtree by specifying the node
##'     number along with \code{nodetype = "include"}. Likewise, if you want
##'     to look at just the background rate - excluding some particular
##'     lineage - just specify \code{nodetype = "exclude"}.
##'
##' @return An object of class \code{bamm-ratematrix} with the following
##'     components:
##'
##'     \item{lambda}{A \code{nsamples} x \code{nslices} matrix of speciation
##'         rates, where \code{nsamples} is the number of posterior samples in
##'         the \code{bammdata} object.}
##'     \item{mu}{A \code{nsamples} x \code{nslices} matrix of extinction
##'         rates.}
##'     \item{beta}{A \code{nsamples} x \code{nslices} matrix of phenotypic
##'         rates (if applicable).}
##'     \item{times}{A vector of timepoints where rates were computed.}
##'     \item{times}{A vector of timepoints where rates were computed (see
##'         Examples).}
##'     \item{type}{Either "diversification" or "trait", depending on the
##'         input data.}
##'
##' @author Dan Rabosky
##'
##' @seealso \code{\link{plotRateThroughTime}}
##'
##' @examples
##' \dontrun{
##' # Plot a rate-through-time curve with 
##' # confidence intervals for the whale dataset:
##'
##' data(whales, events.whales)
##' ed <- getEventData(whales, events.whales)
##'
##' rmat <- getRateThroughTimeMatrix(ed)
##'
##' plot.new()
##' plot.window(xlim=c(0, 36), ylim=c(0, .7))
##'
##' ## Speciation quantiles: plot 90% CIs
##' qq <- apply(rmat$lambda, 2, quantile, c(0.05, 0.5, 0.95))
##'
##' xv <- c(rmat$times, rev(rmat$times))
##' yv <- c(qq[1,], rev(qq[3,]))
##'
##' ## Add the confidence polygon on rate distributions:
##' polygon(xv, yv, col='gray80', border=FALSE)
##'
##' ## Add the median rate line:
##' lines(rmat$times, qq[2,], lwd=3, col='red')
##'
##' ## Add axes
##' axis(1, at=seq(-5, 35, by=5))
##' axis(2, at=seq(-0.2, 1, by=0.2), las=1)
##'
##' ####### Now we will show the actual grid used for rate calculations:
##'
##' plot(whales, show.tip.label=FALSE)
##' axisPhylo()
##'
##' mbt <- max(branching.times(whales))
##' tvec <- mbt - rmat$times;
##' tvec <- rmat$times;
##'
##' for (i in 1:length(tvec)){
##'     lines(x=c(tvec[i], tvec[i]), y=c(0, 90), lwd=0.7, col='gray70')
##' }
##'
##' ## This shows the grid of time slices over the phylogeny}
##' @keywords models
##' @export
getRateThroughTimeMatrix <- function(ephy, start.time=NULL, end.time=NULL, nslices=100, node=NULL, nodetype = 'include') {
	
	if (!inherits(ephy, 'bammdata')) {
		stop("Object ephy must be of class 'bammdata'\n");
	}
	
	if (is.null(node)) {
		nodeset <- c(length(ephy$tip.label) + 1, ephy$edge[,2]);
	} else if (!is.null(node) & nodetype == 'include') {
		nodeset <- unlist(sapply(node, function(x) getDesc(ephy, x)$desc_set))
	} else if (!is.null(node) & nodetype == 'exclude') {
		nodeset <- setdiff( ephy$edge[,2], unlist(sapply(node, function(x) getDesc(ephy, x)$desc_set)));
		nodeset <- c(length(ephy$tip.label) + 1, nodeset);
	} else {
		stop('error in getRateThroughTimeMatrix\n');
	}
	
	if (is.ultrametric(as.phylo.bammdata(ephy))) {
		bt <- branching.times(as.phylo.bammdata(ephy));
	}
	if (!is.ultrametric(as.phylo.bammdata(ephy))) {
		bt <- NU.branching.times(as.phylo.bammdata(ephy));
	}
	maxpossible <- max(bt[as.character(intersect(nodeset, ephy$edge[,1]))]);


	#convert from time before present to node heights
	if (!is.null(start.time)) {
		new.start.time <- max(bt) - start.time;
	}
	if (!is.null(end.time)) {
		new.end.time <- max(bt) - end.time;
	}


	if (is.null(start.time)) {
		new.start.time <- max(bt) - maxpossible;
		start.time <- maxpossible;
	}
	if (is.null(end.time)) {
		new.end.time <- max(bt);
		end.time <- 0;
	}
	
	tvec <- seq(new.start.time, new.end.time, length.out = nslices);
	names(tvec) <- seq(start.time, end.time, length.out = nslices);
	#tol = 1*10^-decimals(ephy$eventBranchSegs[[1]][1,2]);
	tol <- 0.00001
		
	mm <- matrix(NA, nrow = length(ephy$eventBranchSegs), ncol = length(tvec));
	mumat <- matrix(NA, nrow = length(ephy$eventBranchSegs), ncol = length(tvec));

	for (i in 1:nrow(mm)){
		es <- ephy$eventBranchSegs[[i]];
		events <- ephy$eventData[[i]];
						
		if (is.null(node)){ 
			isGoodNode <- rep(TRUE, nrow(es));
		} else {
			# es[,1] are descendant nodes of branches. 
			# If es[,1] is in nodeset, then at least part of the branch is in the 
			# tree subset of interest
			isGoodNode <- es[, 1] %in% nodeset;	
		}
		
		# # plot branch segments that are to be included in time binned mean rate calc.
		# plot(as.phylo.bammdata(ephy), cex = 0.5)
		# nodelabels(node = node, frame = 'circle', cex = 0.4)
		# abline(v = tvec, lty = 3, col = gray(0.5))
		# axisPhylo()
		# axis(1, line = 2)

		for (k in 1:length(tvec)){
			isGoodTime <- goodTime(es, tvec[k], tol = tol);
			
			estemp <- matrix(es[isGoodTime & isGoodNode, ], ncol = 4);
			
			# pp <- get("last_plot.phylo", envir = .PlotPhyloEnv)
			# segments(x0 = estemp[,2], x1 = estemp[,3], y0 = pp$yy[estemp[,1]], y1 = pp$yy[estemp[,1]], lwd = 2, col = 'blue')
			
			tvv <- tvec[k] - events$time[estemp[,4]];
			rates <- exponentialRate(tvv, events$lam1[estemp[,4]], events$lam2[estemp[,4]]);
			mm[i, k] <- mean(rates);
			mumat[i,k] <- mean(events$mu1[estemp[,4]]);	
		}	
	}

	obj <- list();
	if (ephy$type == 'diversification') {
		obj$lambda <- mm;
		obj$mu <- mumat;
	}
	if (ephy$type == 'trait') {
		obj$beta <- mm;
	}
	obj$times <- tvec;
	
	class(obj) <- 'bamm-ratematrix';
	if (ephy$type == 'diversification') {
		obj$type <- 'diversification';
	} else {
		obj$type <- 'trait';	
	}
	return(obj);
}



goodTime <- function (vec, val, tol) {
	(vec[,2] - val <= tol) & (val - vec[,3] <= tol);
}

