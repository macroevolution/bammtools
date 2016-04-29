##' @title Summary of rate shift results from \code{BAMM} analysis
##'
##' @description Summarizes the posterior distribution on the number of
##'     shifts.
##'
##' @param object An object of class \code{bammdata}.
##' @param display An integer for the number of rows of the posterior to
##'     display.
##' @param print Print summary of shift distribution in console window?
##' @param \dots Additional arguments (currently unused).
##'
##' @details Prints to console the number of posterior samples and the
##'     posterior distribution on the number of shifts, which is just the
##'     fraction of samples in the posterior having 0, 1, 2,...n shifts.
##'
##' @return Returns (invisibly) a dataframe with 2 components:
##'     \item{shifts}{The number of shifts.}
##'     \item{prob}{The corresponding posterior probability of a model with a
##'     given number of rate shifts.}
##'
##' @author Mike Grundler, Dan Rabosky
##'
##' @references \url{http://bamm-project.org}
##'
##' @examples
##' data(whales, events.whales)
##' ephy <- getEventData(whales, events.whales, nsamples=100)
##' summary(ephy)
##' @keywords models
##' @export
summary.bammdata = function(object, display=10, print=T, ...) {

	fev <- sapply(object$eventData, nrow);
	xx <- table(fev - 1);
	shifts <- as.numeric(names(xx));
	xx <- as.numeric(xx);
	df <- data.frame(shifts = shifts, prob = signif(xx / sum(xx), digits= 3));

	if (print){
		
		cat("\nAnalyzed", length(object$eventData), "posterior samples\n");
	 
		cat("Shift posterior distribution:\n\n");
		disp <- data.frame(cbind(df$shifts),signif(df$prob,2));
		if (nrow(disp) <= display) {	
			write.table(format(disp,justify="left",width=10), col.names=FALSE, row.names=FALSE, quote=FALSE);	
		} else {
			if(nrow(disp)-display == 1) {
				write.table(format(disp[1:display,],justify="left",width=10), col.names=FALSE, row.names=FALSE, quote=FALSE);
				cat("... omitted 1 row\n\n");
			}else {
				wr <- which(disp[,2] == max(disp[,2]))[1];
				index <- c(max(1,floor(wr-display/2)), wr, min(ceiling(wr+display/2),nrow(disp)));
				if (index[1] > 1) {
					cat("... omitted",index[1]-1,"rows\n");
				}
				write.table(format(disp[index[1]:index[3],],justify="left",width=10), col.names=FALSE, row.names=FALSE, quote=FALSE);
				cat("... omitted", nrow(disp)-index[3]-1,"rows\n\n");
			}
		}
		cat("\nCompute credible set of shift configurations for more information:\n");
		cat("\tSee ?credibleShiftSet and ?getBestShiftConfiguration\n");		
		
	}
	
	invisible(df);
}
