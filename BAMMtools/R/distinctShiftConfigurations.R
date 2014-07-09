


# Drop all nodes from event data with marginal probs < 0.05
# test is unique
#  if so, add to list
# Returns:
#	$marg.probs = marginal probs for nodes
#	$bayesfactors = branch-specific BF for shift
#	$shifts = unique shift sets
#	$samplesets = list of sample indices that reduce to each of the unique shift sets
#	$frequency = vector of frequencies of each shift configuration
#	$bfthreshold = bayes factor threshold for shifts
#	
#	Results are sorted by frequency. 
#	$frequency[1] gives the most common shift config sampled
#	$shifts[[1]] gives the corresponding node indices for that configuration
#	$samplesets[[1]] gives the indices of samples with this configuration



distinctShiftConfigurations <- function(ephy, prior, BFcriterion, ... ) {
	
	if (hasArg("threshold")){
		cat("Argument < threshold > has been deprecated. It is \nreplaced");
		cat(" by the argument < BFcriterion >, \nwhich uses an explicit Bayes factor")
		cat(" criterion to identify core shifts.\n Please see help on this function")
		cat(" ( ?distinctShiftConfigurations )\n\n");
		cat("Apologies for the change, but the new way is much better...\n")
		stop();
		
	}
	
	if (class(prior) != 'branchprior'){
		stop("object prior not of class branchprior");
	}
	
	bf <- bayesFactorBranches(ephy, prior);
	
	mm <- marginalShiftProbsTree(ephy);

	goodnodes <- bf$edge[,2][bf$edge.length >= BFcriterion];

	xlist <- list();
	for (i in 1:length(ephy$eventData)) {
		xlist[[i]] <- intersect(goodnodes, ephy$eventData[[i]]$node);
	}

	ulist <- list();
	treesets <- list();
	
	ulist[[1]] <- xlist[[1]];
	treesets[[1]] <- 1;
	
	for (i in 2:length(xlist)) {
		lx <- length(ulist);
		#cat(lx, '\n')
		for (k in 1:lx) {
			if (areShiftSetsEqual(ulist[[k]], xlist[[i]])){
				treesets[[k]] <- c(treesets[[k]], i);
				break;	
			} else {
				if (k == length(ulist)){
					xlen <- length(ulist);
					ulist[[xlen + 1]] <- xlist[[i]];
					treesets[[xlen + 1]] <- i;
				}
			}
		}
	}
	
	freqs <- unlist(lapply(treesets, length));
	freqs <- freqs / sum(freqs);
	
	ord <- order(freqs, decreasing=TRUE);
	
	obj <- list();
	obj$marg.probs <- mm$edge.length;  
	names(obj$marg.probs) <- mm$edge[,2]; 
	obj$bayesfactors <- bf$edge.length;
	names(obj$bayesfactors) <- bf$edge[,2];
	obj$shifts <- ulist[ord]; 
	obj$samplesets <- treesets[ord];
	obj$frequency <- freqs[ord];
	obj$coreshifts <- goodnodes;
	obj$BFcriterion <- BFcriterion;

	class(obj) <- 'bammshifts';
	
	return(obj);
}






