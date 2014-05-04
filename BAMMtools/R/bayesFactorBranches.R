
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
	tree$edge.length <- mprobs$edge.length / prior$edge.length;
	return(tree);
}





