##' @title Add \code{BAMM}-inferred rate shifts to a phylogeny plot
##'
##' @description Adds symbols to a plotted tree to mark the location(s) where
##'     there is a shift in the macroevolutionary dynamics of diversification
##'     or trait evolution.
##'
##' @param ephy An object of class \code{bammdata}.
##' @param method A character string indicating the method used in plotting.
##'     Must be "polar" or "phylogram".
##' @param index An integer indicating which posterior sample to use for
##'     adding shifts to the plotted tree.
##' @param cex A numeric indicating the character expansion ("size") of the
##'     plotted points.
##' @param pch An integer indicating the choice of plotting symbol.
##' @param col An integer or character string indicating the border color of
##'     the plotting symbol.
##' @param bg An integer or character string indicating the background color
##'     of the plotting symbol.
##' @param msp If not \code{NULL}, an object of class \code{phylo} where each
##'     branch length is equal to the marginal probability of a shift
##'     occurring on that branch. Plotted points corresponding to shifts will
##'     be sized by these probabilities.  	
##' @param shiftnodes An optional vector of node numbers indicating the
##'     locations of shifts to plot. 	
##' @param par.reset A logical indicating whether to reset the graphical
##'     parameters before exiting.
##'
##' @details Any given sample from the posterior distribution sampled using
##'     \code{BAMM} contains a potentially unique configuration of rate shifts
##'     and associated parameters. There is no single "best" rate shift, but
##'     rather a set of shift configurations (and associated parameters) -
##'     along with their relative probabilities - sampled with MCMC. This
##'     function enables the user to plot the locations of shifts sampled with
##'     \code{BAMM} for a given sample from the posterior. 
##'
##'     If the \code{bammdata} object contains just a single sample, these
##'     shifts will be plotted regardless of the value of \code{index}.
##'
##' @note If a \code{shiftnodes} argument is passed care should be taken to
##'     ensure that the nodes are in the same order as in the event data for
##'     the sample index.
##'
##' @author Mike Grundler
##'
##' @seealso \code{\link{getShiftNodesFromIndex}}, \code{\link{plot.bammdata}}
##'
##' @examples
##' data(whales, events.whales)
##' ed <- getEventData(whales, events.whales, burnin=0.25, nsamples = 500)
##' 
##' # adding shifts to tree for specific posterior samples
##' plot(ed, method="polar")
##' addBAMMshifts(ed, index=5, "polar")
##' 
##' # multi-panel plotting and adding shifts
##' par(mfrow=c(2,3),mar=c(5,1,1,1))
##' samples = sample(1:length(ed$eventData), 6)
##' for (i in 1:6) {
##'   sed <- subsetEventData(ed, samples[i])
##'   plot(sed, par.reset=FALSE)
##'   addBAMMshifts(sed,index=1,method="phylogram",par.reset=FALSE)	
##' }
##' @keywords graphics
##' @export
addBAMMshifts = function(ephy, index = 1, method = 'phylogram', cex=1, pch=21, col=1, bg=2, msp = NULL, shiftnodes = NULL, par.reset=TRUE) {
	if (!'bammdata' %in% class(ephy)) stop("Object ephy must be of class bammdata");
	lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv);
	
	if (par.reset){
		op <- par(no.readonly = TRUE);
		par(lastPP$pp);		
	}

	if (length(ephy$eventData) == 1){
		index <- 1;
	}
	
	if (is.null(shiftnodes))
		shiftnodes <- getShiftNodesFromIndex(ephy, index)
	isShift <- ephy$eventData[[index]]$node %in% shiftnodes;
	times <- ephy$eventData[[index]]$time[isShift];	
	if (!is.null(msp)) {
		cex <- 0.75 + 5 * msp$edge.length[msp$edge[,2] %in% shiftnodes];
	}
	
	if (method == 'phylogram') {
		###  obsolete b/c plot.bammdata no longer scales each axis to a max of 1. now behaves like plot.phylo
		# if (max(lastPP$xx) <= 1) {
		# 	XX <- times/max(branching.times(as.phylo.bammdata(ephy)));
		# } else {
		# 	XX <- times;
		# }
		XX <- times;
		YY <- lastPP$yy[shiftnodes];
	} else if (method == 'polar') {
		rb <- lastPP$rb;
		XX <- (rb+times/max(branching.times(as.phylo.bammdata(ephy)))) * cos(lastPP$theta[shiftnodes]);
		YY <- (rb+times/max(branching.times(as.phylo.bammdata(ephy)))) * sin(lastPP$theta[shiftnodes]);		
	}	
	points(XX,YY,pch=pch,cex=cex,col=col,bg=bg);
	if (par.reset) {
		par(op);		
	}
}
