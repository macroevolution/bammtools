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
	
	if (!'bammdata' %in% class(ephy)) {
		stop("Object ephy must be of class 'bammdata'\n");
	}
	
	if (is.null(node)) {
		nodeset <- c(length(ephy$tip.label)+1, ephy$edge[,2]);
	} else if (!is.null(node) & nodetype == 'include') {
		nodeset <- unlist(sapply(node, function(x) getDesc(ephy, x)$desc_set))
	} else if (!is.null(node) & nodetype == 'exclude') {
		nodeset <- setdiff( ephy$edge[,2], unlist(sapply(node, function(x) getDesc(ephy, x)$desc_set)));
		nodeset <- c(length(ephy$tip.label)+1, nodeset);
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
		start.time <- max(bt) - start.time;
	}
	if (!is.null(end.time)) {
		end.time <- max(bt) - end.time;
	}


	if (is.null(start.time)) {
		start.time <- max(bt) - maxpossible;
	}
	if (is.null(end.time)) {
		end.time <- max(bt);
	}
	
 	#remove node from nodeset to prevent branch leading up to node from being included
	# This is only an issue if start.time occurs before node.
	if (!is.null(node)) {
		if (nodetype == 'include') {
			nodeset <- nodeset[nodeset != node];
		}
	}

	
	tvec <- seq(start.time, end.time, length.out= nslices);
	#tol = 1*10^-decimals(ephy$eventBranchSegs[[1]][1,2]);
	tol <- 0.00001
	
	goodTime <- function (vec, val, tol) {
		(vec[,2] - val <= tol) & (val - vec[,3] <= tol);
	}
	
	getRates <- function(time, es, events, isGoodNode) {
		isGoodTime <- goodTime(es, time, tol=tol);
		if (!(is.null(node))) { # only enter this if not the root. otherwise, only have to set once per i.
			isGoodNode <- es[,1] %in% nodeset;	
		}
		estemp <- es[isGoodTime & isGoodNode, ];
		if (is.vector(estemp)) {
			index <- estemp[4];
		} else {
			index <- estemp[,4];
		}
		tvv <- time - events$time[index];
		rates <- exponentialRate(tvv, events$lam1[index], events$lam2[index]);
		return(list(rates,index));
	}

	bySample <- function(counter, ephy) {
		es <- ephy$eventBranchSegs[[counter]];
		events <- ephy$eventData[[counter]];
		isGoodNode <- rep(TRUE, nrow(es));
		ret <- lapply(tvec, function(x) getRates(time = x, es, events, isGoodNode));
		mmRow <- unlist(lapply(ret, function(x) mean(x[[1]])));
		mumatRow <- unlist(lapply(ret, function(x) mean(events$mu1[x[[2]]])));
		return(list(mmRow,mumatRow));
	}

	
	ret <- lapply(1:length(ephy$eventBranchSegs), function(y) bySample(y, ephy));
	mm <- lapply(ret, function(x) x[[1]]);
	mm <- do.call(rbind, mm);
	mumat <- lapply(ret, function(x) x[[2]]);
	mumat <- do.call(rbind, mumat);
		
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
	if (ephy$type=='diversification') {
		obj$type = 'diversification';
	} else {
		obj$type = 'trait';	
	}
	return(obj);
}




# # # non-apply version for debugging

# getRateThroughTimeMatrix <- function(ephy, start.time=NULL, end.time=NULL, nslices=100, node=NULL, nodetype = 'include') {
	
	# if (!'bammdata' %in% class(ephy)) {
		# stop("Object ephy must be of class 'bammdata'\n");
	# }
	
	# if (is.null(node)) {
		# nodeset <- c(length(ephy$tip.label)+1, ephy$edge[,2]);
	# } else if (!is.null(node) & nodetype == 'include') {
		# nodeset <- unlist(sapply(node, function(x) getDesc(ephy, x)$desc_set))
	# } else if (!is.null(node) & nodetype == 'exclude') {
		# nodeset <- setdiff( ephy$edge[,2], unlist(sapply(node, function(x) getDesc(ephy, x)$desc_set)));
		# nodeset <- c(length(ephy$tip.label)+1, nodeset);
	# } else {
		# stop('error in getRateThroughTimeMatrix\n');
	# }
	
	# if (is.ultrametric(as.phylo.bammdata(ephy))) {
		# bt <- branching.times(as.phylo.bammdata(ephy));
	# }
	# if (!is.ultrametric(as.phylo.bammdata(ephy))) {
		# bt <- NU.branching.times(as.phylo.bammdata(ephy));
	# }
	# maxpossible <- max(bt[as.character(intersect(nodeset, ephy$edge[,1]))]);

	# #convert from time before present to node heights
	# if (!is.null(start.time)) {
		# start.time <- max(bt) - start.time;
	# }
	# if (!is.null(end.time)) {
		# end.time <- max(bt) - end.time;
	# }

	# if (is.null(start.time)) {
		# start.time <- max(bt) - maxpossible;
	# }
	# if (is.null(end.time)) {
		# end.time <- max(bt);
	# }
	
	# tvec <- seq(start.time, end.time, length.out= nslices);
	# #tol = 1*10^-decimals(ephy$eventBranchSegs[[1]][1,2]);
	# tol <- 0.00001
	
	# goodTime <- function (vec, val, tol) {
		# (vec[,2] - val <= tol) & (val - vec[,3] <= tol);
	# }	
	
	# #remove node from nodeset to prevent branch leading up to node from being included
	# # This is only an issue if start.time occurs before node.
	# if (!is.null(node)) {
		# if (nodetype == 'include') {
			# nodeset <- nodeset[nodeset != node];
		# }
	# }

	
	# mm <- matrix(NA, nrow=length(ephy$eventBranchSegs), ncol=length(tvec));
	# mumat <- matrix(NA, nrow=length(ephy$eventBranchSegs), ncol=length(tvec));

	# for (i in 1:nrow(mm)){
		# es <- ephy$eventBranchSegs[[i]];
		# events <- ephy$eventData[[i]];
				
		# for (k in 1:length(tvec)){
			# isGoodTime <- goodTime(es, tvec[k], tol=tol);
			
			# if (is.null(node)){ 
				# isGoodNode <- rep(TRUE, nrow(es));
			# } else {
				# isGoodNode <- es[,1] %in% nodeset;	
			# }
			# estemp <- es[isGoodTime & isGoodNode, ];
			# tvv <- tvec[k] - events$time[estemp[,4]];
			# rates <- exponentialRate(tvv, events$lam1[estemp[,4]], events$lam2[estemp[,4]]);
			# mm[i, k] <- mean(rates);
			# mumat[i,k] <- mean(events$mu1[estemp[,4]]);	
		# }	
	# }
			
	# obj <- list();
	# if (ephy$type == 'diversification') {
		# obj$lambda <- mm;
		# obj$mu <- mumat;
	# }
	# if (ephy$type == 'trait') {
		# obj$beta <- mm;
	# }
	# obj$times <- tvec;
	
	# class(obj) <- 'bamm-ratematrix';
	# if (ephy$type=='diversification') {
		# obj$type = 'diversification';
	# } else {
		# obj$type = 'trait';	
	# }
	# return(obj);
# }
