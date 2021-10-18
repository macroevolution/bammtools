##' @title Compute species-specific rate through time trajectories
##'
##' @description Computes the mean of the marginal posterior density of
##'     speciation/extinction or phenotypic rates for equally spaced points
##'     along the root to tip path for each species.
##'
##' @param ephy An object of class \code{bammdata}.
##' @param nslices An integer number of time slices. This determines the
##'     number of equally spaced points in time at which rates are computed
##'     for each species.
##' @param index An integer or vector of mode integer indicating which
##'     posterior samples to use in the calculation. If \code{NULL} (default)
##'     all samples are used.
##' @param spex A character string. "s" (default) calculates speciation rates;
##'     "e" calculates extinction rates; "netdiv" calculates diversification
##'     rates. Ignored if \code{ephy$type = "trait"}.
##'
##' @return A list with two components:
##'     \itemize{
##'         \item{times} {A vector of time points where rates were
##'             calculated.}
##'         \item{rates} {A species X times matrix of rate through time
##'             trajectories.}
##'     }
##'
##' @author Mike Grundler
##'
##' @seealso \code{\link{getRateThroughTimeMatrix}}
##'
##' @references \url{http://bamm-project.org/}
##'
##' @examples
##' data(whales, events.whales)
##' ed <- getEventData(whales, events.whales, burnin=0.25, nsamples=500)
##' ratemat <- speciesByRatesMatrix(ed, nslices = 100)
##' 
##' dolphins <- extract.clade(whales, 140)$tip.label
##' plot.new()
##' plot.window(xlim = c(0, 35), ylim = c(0, 0.8))
##' for (i in 1:nrow(ratemat$rates)) {
##' 	if (whales$tip.label[i] %in% dolphins) {
##' 		lines(ratemat$times, ratemat$rates[i,], lwd = 2, col = 4)	
##' 	} else {
##' 		lines(ratemat$times, ratemat$rates[i,], lwd = 2, col = 8)
##' 	}
##' }
##' axis(1, seq(-5, 35, 5))
##' axis(2, seq(-0.2, 0.8, 0.2), las = 1)
##' mtext("Time since root", 1, line = 2.5)
##' mtext("Speciation rate", 2, line = 2.5)
##'
##' @keywords models
##' @export

speciesByRatesMatrix = function(ephy, nslices, index = NULL, spex = "s") {
	if (!spex %in% c('s', 'e', 'netdiv')) {
		stop("arg spex must be 's', 'e' or 'netdiv'.")
	}
	seq.nod <- .Call("seq_root2tip", ephy$edge, length(ephy$tip.label), ephy$Nnode, PACKAGE = "BAMMtools");
	if (nslices <= 100) {
		tvec <- (seq(0, 1, 0.01)+0.005) * max(ephy$end);
		tvec <- tvec[seq.int(1,length(tvec),length.out=nslices+1)];
		ephy <- dtRates(ephy, 0.01, index, tmat = TRUE);
	}
	else if (nslices > 100 && nslices <= 500) {
		tvec <- (seq(0, 1, 0.002)+0.001) * max(ephy$end);
		tvec <- tvec[seq.int(1,length(tvec),length.out=nslices+1)];
		ephy <- dtRates(ephy, 0.002, index, tmat = TRUE);
	}
	else if (nslices > 500 && nslices <= 1000) {
		tvec <- (seq(0, 1, 0.001)+0.0005) * max(ephy$end);
		tvec <- tvec[seq.int(1,length(tvec),length.out=nslices+1)];
		ephy <- dtRates(ephy, 0.001, index, tmat = TRUE);
	}
	else {
		stop("Max slices (1000) exceeded.  Choose a smaller number of slices");
	}
	ret <- lapply(seq.nod, function(x) {
		path = which(ephy$dtrates$tmat[,1] %in% x);
		ids = sapply(tvec[-length(tvec)], function(y) which(ephy$dtrates$tmat[path,2] <= y & ephy$dtrates$tmat[path,3] > y));
		if (is.list(ids))
			ids = unlist(ids);
		if (ephy$type == "trait") {
			rts = ephy$dtrates$rates[path][ids];
		}
		else {
			if (tolower(spex) == "s") {
				rts = ephy$dtrates$rates[[1]][path][ids];
			}
			else if (tolower(spex) == "e") {
				rts = ephy$dtrates$rates[[2]][path][ids];
			}
			else if (tolower(spex) == "netdiv") {
				rts = ephy$dtrates$rates[[1]][path][ids] - ephy$dtrates$rates[[2]][path][ids];
			}
		}
		if (length(rts) < (length(tvec)-1))
			rts = c(rts, rep(NA, length(tvec)-1-length(rts)));
		rts;
	});
	ret <- do.call(rbind, ret);
	rownames(ret) <- ephy$tip.label;
	return(list(times = tvec[-length(tvec)],rates = ret));	
}
