#############################################################
#
#	getCladeRates(....)
#
#	mean clade-specific rates 
#		average of all branch-rates, but weighted by branch length
#		node.type: will compute rates only for clade descended from specified node with 'include'
#					will compute for all branches excluding a given clade, nodetype = 'exclude'
#		

##' @title Compute clade-specific mean rates
##'
##' @description Computes marginal clade-specific rates of speciation,
##'     extinction, or (if relevant) trait evolution from \code{BAMM} output.
##'
##' @param ephy An object of class \code{bammdata}.
##' @param node If computing rates for a specific portion of tree, the node
##'     subtending the relevant subtree. If multiple nodes are supplied, then
##'     the equivalent number of \code{nodetype} must be supplied.
##' @param nodetype Either "include" or "exclude". If \code{nodetype = 
##'     "include"}, the rates returned by the function will be for the subtree
##'     defined by \code{node}. If \code{nodetype = "exclude"}, will compute
##'     mean rates for the tree after excluding the subtree defined by
##'     \code{node}. If multiple nodes are specified, there must be a
##'     \code{nodetype} for each node.
##' @param verbose Logical. If \code{TRUE}, will print the sample index as
##'     mean rates are computed for each sample from posterior. Potentially
##'     useful for extremely large trees.
##'
##' @details Computes the time-weighted mean evolutionary rate for a given
##'     clade. Conversely, one can compute the rate for a given phylogeny
##'     while excluding a clade; this operation will give the "background"
##'     rate. It is important to understand several aspects of these mean
##'     rates. First, rates in the \code{BAMM} framework are not constant
##'     through time. Hence, the function computes the mean time-integrated
##'     rates across the subtree. Operationally, this is done by integrating
##'     the speciation rate with respect to time along each branch in the
##'     subtree. These time-integrated rates are then summed, and the sum
##'     is divided by the total sum of branch lengths for the subtree. 
##'
##'     The function computes a rate for each sample in the posterior, and
##'     returns a list of rate vectors. Each rate in the corresponding vector
##'     is a mean rate for a particular sample from the posterior. Hence, you
##'     can think of the return value for this function as an estimate of the
##'     marginal distribution of rates for the focal clade. You can compute
##'     means, medians, quantiles, etc from these vectors.
##'
##' @return A list with the following components:
##'     \itemize{
##'         \item{lambda} {A vector of speciation rates (if applicable), with
##'             the i'th rate corresponding to the mean rate from the i'th
##'             sample in the posterior.}
##'         \item{mu} {A vector of extinction rates (if applicable), with the
##'             i'th rate corresponding to the mean rate from the i'th sample
##'             in the posterior.}
##'         \item{beta} {A vector of phenotypic rates (if applicable), with
##'             the i'th rate corresponding to the mean rate from the i'th 
##'             sample in the posterior.}
##'     }
##'
##' @author Dan Rabosky
##'
##' @references \url{http://bamm-project.org/}
##'
##' @examples
##' data(events.whales, whales)
##' ed <- getEventData(whales, events.whales, nsamples=500)
##' all_rates <- getCladeRates(ed)
##' 
##' mean(all_rates$lambda)
##' mean(all_rates$mu)
##' # joint density of mean speciation and extinction rates:
##' plot(all_rates$mu ~ all_rates$lambda)
##' 
##' # clade specific rates: here for Dolphin subtree:
##' dol_rates <- getCladeRates(ed, node=140)
##' mean(dol_rates$lambda)
##' mean(dol_rates$mu)
##' 
##' # defining multiple nodes
##' mean(getCladeRates(ed, node=c(132, 140),
##'      nodetype=c('include','exclude'))$lambda)
##' @keywords models
##' @export
getCladeRates <- function(ephy, node = NULL, nodetype='include', verbose=FALSE) {
	
	if (!inherits(ephy, 'bammdata')) {
		stop("Object ephy must be of class bammdata\n");
	}	
	
	if (is.null(node)) {
		nodeset <- ephy$edge[,2];
	} else if (!is.null(node[1]) & nodetype[1] == 'include' & length(node) == 1) {
		nodeset <- getDesc(ephy, node)$desc_set;
	} else if (!is.null(node[1]) & nodetype[1] == 'exclude' & length(node) == 1) {
		nodeset <- setdiff(ephy$edge[,2],  getDesc(ephy, node)$desc_set);
	} else if (!is.null(node[1]) & length(nodetype) == length(node) & length(node) > 1) {
		nodesets <- lapply(node, function(x) getDesc(ephy, x)$desc_set);
		Drop <- which(nodetype == 'exclude');
		nodeset_toRemove <- unique(unlist(lapply(nodesets[Drop], function(x) x)));
		Keep <- which(nodetype == 'include');
		nodeset_toKeep <- unique(unlist(nodesets[Keep]));
		nodeset <- setdiff(nodeset_toKeep, nodeset_toRemove);
		if (length(nodeset) == 0) {
			stop('Error: the combination of nodes and nodetypes has resulted in no remaining nodes!')
		}
	} else {
		stop('Error: Please make sure you have specified a nodetype for every node');
	}
	
	lambda_vector <- numeric(length(ephy$eventBranchSegs));
	mu_vector <- numeric(length(ephy$eventBranchSegs));
 	
 	weights <- 'branchlengths'
	
	for (i in 1:length(ephy$eventBranchSegs)) {
		if (verbose) {
			cat('Processing sample', i, '\n');
		}
		esegs <- ephy$eventBranchSegs[[i]];
		esegs <- esegs[esegs[,1] %in% nodeset, ];
	
       	if (is.null(nrow(esegs))){
       		esegs <- t(as.matrix(esegs))
       	}	
		
		events <- ephy$eventData[[i]];
		events <- events[order(events$index), ];			
		
		# relative start time of each seg, to event:
		relsegmentstart <- esegs[,2] - ephy$eventData[[i]]$time[esegs[,4]];
		relsegmentend <- esegs[,3] - ephy$eventData[[i]]$time[esegs[,4]];
		lam1 <- ephy$eventData[[i]]$lam1[esegs[,4]];
		lam2 <-  ephy$eventData[[i]]$lam2[esegs[,4]];
		mu1 <-  ephy$eventData[[i]]$mu1[esegs[,4]];
		mu2 <-  ephy$eventData[[i]]$mu2[esegs[,4]];
 		
 		seglengths <- esegs[,3] - esegs[,2];	
		wts <- seglengths / sum(seglengths);
		lamseg <- timeIntegratedBranchRate(relsegmentstart, relsegmentend, lam1, lam2) / seglengths;
		museg <- timeIntegratedBranchRate(relsegmentstart, relsegmentend, mu1, mu2) / seglengths;
	
		lambda_vector[i] <- sum(lamseg * wts);
		mu_vector[i] <- sum(museg  * wts);
	}
	
	if (ephy$type == 'diversification') {
		return(list(lambda = lambda_vector, mu = mu_vector));
	}
	if (ephy$type == 'trait') {
		return(list(beta = lambda_vector));
	}
}
