
# returns a phylogenetic tree where 
# branch lengths are equal to the Bayes factor
# evidence in favor of a rate shift.

bayesFactorBranches <- function(ephy, priordata){
	
	mprobs <- marginalShiftProbsTree(ephy);
	if (class(priordata) == 'character'){
		prior <- getBranchShiftPriors(prior);			
	} else if (class(priordata) == 'branchprior'){
		prior <- priordata;
	}else{
		stop('object priordata of wrong class\n');
	}

	tree <- mprobs;
	
	post_shift <- mprobs$edge.length;
	prior_shift <- prior$edge.length;
	
	post_noshift <- 1 - mprobs$edge.length;
	prior_noshift <- 1 - prior$edge.length;
	
	bf <- (post_shift / post_noshift) * (prior_noshift / prior_shift);

	tree$edge.length <- bf;
	return(tree);
}





