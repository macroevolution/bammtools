##' @title Get the best (sampled) rate shift configuration from a \code{BAMM}
##'     analysis
##'
##' @description Get the rate shift configuration with the maximum a
##'     posteriori probability, e.g., the shift configuration that was sampled
##'     most frequently with \code{BAMM}.
##'
##' @param x Either a \code{bammdata} object or a \code{credibleshiftset}
##'     object.
##' @param expectedNumberOfShifts The expected number of shifts under the
##'     prior.
##' @param threshold The marginal posterior-to-prior odds ratio used as a
##'     cutoff for distinguishing between "core" and "non-core" rate shifts.
##'
##' @details This function estimates the rate shift configuration with the
##'     highest maximum a posteriori (MAP) probability. It returns a
##'     \code{bammdata} object with a single sample. This can be plotted with
##'     \code{\link{plot.bammdata}}, and individual rate shifts can then
##'     be added with \code{\link{addBAMMshifts}}. 
##'
##'     The parameters of this object are averaged over all samples in the
##'     posterior that were assignable to the MAP shift configuration. All
##'     non-core shifts have been excluded, such that the only shift
##'     information contained in the object is from the "significant" rate
##'     shifts, as determined by the relevant marginal posterior-to-prior odds
##'     ratio \code{threshold}.
##'
##'     You can extract the same information from the credible set of shift
##'     configurations. See \code{\link{credibleShiftSet}} for more
##'     information.
##'
##' @return A class \code{bammdata} object with a single sample, corresponding
##'     to the diversification rate shift configuration with the maximum a
##'     posteriori probability. See \code{\link{getEventData}} for details.
##'
##' @author Dan Rabosky
##'
##' @seealso \link{getEventData},  \link{credibleShiftSet},
##'     \link{plot.credibleshiftset}, \link{plot.bammdata}
##'
##' @examples
##' data(whales, events.whales)
##' ed <- getEventData(whales, events.whales, burnin=0.1, nsamples=500)
##' 
##' # Get prior distribution on shifts-per-branch:
##' bp <- getBranchShiftPriors(whales, expectedNumberOfShifts = 1)
##' 
##' # Pass the event data object in to the function:
##' best <- getBestShiftConfiguration(ed, expectedNumberOfShifts = 1,
##'                                   threshold = 5)
##' plot(best, lwd=2)
##' addBAMMshifts(best, cex=2)
##' 
##' # Now we can also work with the credible shift set:
##' css <- credibleShiftSet(ed, expectedNumberOfShifts = 1, threshold = 5)
##' 
##' summary(css)
##' 
##' # examine model-averaged shifts from MAP configuration-
##' # This gives us parameters, times, and associated nodes
##' #   of each evolutionary rate regime (note that one of
##' #   them corresponds to the root)
##' css$eventData[[1]];
##' 
##' # Get bammdata representation of MAP configuration:
##' best <- getBestShiftConfiguration(css, expectedNumberOfShifts = 1,
##'                                   threshold = 5)
##' 
##' plot(best)
##' addBAMMshifts(best)
##' @keywords models
##' @export
getBestShiftConfiguration <- function(x, expectedNumberOfShifts , threshold = 5){
	
	if (inherits(x, 'bammdata')) {
		x <- credibleShiftSet(x, expectedNumberOfShifts, threshold, set.limit = 0.95);	
	} else if (inherits(x, 'credibleshiftset')) {

	} else {
		stop("Argument x must be of class bammdata or credibleshiftset\n");
	}
	
	class(x) <- 'bammdata';	
	subb <- subsetEventData(x, index = x$indices[[1]]);
	
	# Drop all non-core shifts after adding root:
	coreshifts <- c((length(x$tip.label) + 1), x$coreshifts);
	coreshifts <- intersect(subb$eventData[[1]]$node, coreshifts); 
	 
		
 	for (i in 1:length(subb$eventData)) {
 		#subb$eventData[[i]] <- subb$eventData[[i]][subb$eventData[[i]]$node %in% coreshifts,];
 		if (i == 1) {
 			ff <- subb$eventData[[i]];
 		}
 		ff <- rbind(ff, subb$eventData[[i]]);
 	}

	xn <- numeric(length(coreshifts));
	xc <- character(length(coreshifts));
	
	if (x$type == 'diversification') {
		dff <- data.frame(generation = xn, leftchild=xc, rightchild=xc, abstime=xn, lambdainit=xn, lambdashift=xn, muinit = xn, mushift = xn, stringsAsFactors=F);	
		for (i in 1:length(coreshifts)) {
			if (coreshifts[i] <= length(x$tip.label)) {
			# Node is terminal:	
				dset <- c(x$tip.label[coreshifts[i]], NA)
				
			}else{
				# node is internal.
				tmp <- extract.clade(as.phylo(x), node= coreshifts[i]);
				dset <- tmp$tip.label[c(1, length(tmp$tip.label))];	
			}
			
			tmp2 <- ff[ff$node == coreshifts[i], ];
			
			dff$leftchild[i] <- dset[1];
			dff$rightchild[i] <- dset[2];
			dff$abstime[i] <- mean(tmp2$time);
			dff$lambdainit[i] <- mean(tmp2$lam1);
			dff$lambdashift[i] <- mean(tmp2$lam2);
			dff$muinit[i] <- mean(tmp2$mu1);
			dff$mushift[i] <- mean(tmp2$mu2);
		}	
		best_ed <- getEventData(as.phylo(x), eventdata=dff);	
	} else if (x$type == 'trait') {
		dff <- data.frame(generation = xn, leftchild=xc, rightchild=xc, abstime=xn, betainit=xn, betashift=xn, stringsAsFactors=F);					
		for (i in 1:length(coreshifts)) {
			if (coreshifts[i] <= length(x$tip.label)) {
			# Node is terminal:	
				dset <- c(x$tip.label[coreshifts[i]], NA)
				
			} else {
				# node is internal.
				tmp <- extract.clade(as.phylo(x), node= coreshifts[i]);
				dset <- tmp$tip.label[c(1, length(tmp$tip.label))];	
			}
			
			tmp2 <- ff[ff$node == coreshifts[i], ];
			
			dff$leftchild[i] <- dset[1];
			dff$rightchild[i] <- dset[2];
			dff$abstime[i] <- mean(tmp2$time);
			dff$betainit[i] <- mean(tmp2$lam1);
			dff$betashift[i] <- mean(tmp2$lam2);

		}			
		best_ed <- getEventData(as.phylo(x), eventdata=dff, type = 'trait');	
	} else {
		stop("error in getBestShiftConfiguration; invalid type");
	}
	return(best_ed);
}

