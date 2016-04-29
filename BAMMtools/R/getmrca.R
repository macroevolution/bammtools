##' @title Find most recent common ancestors
##'
##' @description Calculates the most recent common ancestor for each pair of
##'     tips. Used internally by \code{\link{getEventData}}.
##'
##' @param phy An object of class \code{phylo}.
##' @param t1 A vector of mode integer or character corresponding to tips in
##'     \code{phy}.
##' @param t2 A vector of mode integer or character corresponding to tips in
##'     \code{phy}.
##'
##' @details Finds the most recent common ancestor for each pair of tips where
##'     pairs are defined as (\code{t1}[1], \code{t2}[1]), (\code{t1}[2],
##'     \code{t2}[2]), ... , (\code{t1}[i], \code{t2}[i]), ... ,(\code{t1}[n],
##'     \code{t2}[n]).
##'
##' @return A vector of node numbers of the common ancestor for each pair of
##'     tips.
##'
##' @author Mike Grundler
##'
##' @seealso \code{\link{subtreeBAMM}}
##' @keywords manip
##' @export
getmrca <- function(phy,t1,t2)
{
	if (mode(t1) == "character") {
		t1 <- match(t1, phy$tip.label);
	}
	if (mode(t2) == "character") {
		t2 <- match(t2, phy$tip.label);
	}
	ne <- as.integer(dim(phy$edge)[1]);
	npair <- as.integer(length(t1));
	anc <- as.integer(phy$edge[,1]);
	desc <- as.integer(phy$edge[,2]);
	root <- as.integer(length(phy$tip.label) + 1);
	
	.C('fetchmrca',anc,desc,root,ne,npair,as.integer(t1),as.integer(t2),integer(npair))[[8]];
}
