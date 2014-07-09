
# Feb 28 2014
credibleShiftSet <- function(ephy, prior, BFcriterion = 5, set.limit=0.95, ...){
	
	if (hasArg("threshold")){
		cat("Argument < threshold > has been deprecated. It is \nreplaced");
		cat(" by the argument < BFcriterion >, \nwhich uses an explicit Bayes factor")
		cat(" criterion to identify core shifts.\n Please see help on this function")
		cat(" ( ?credibleShiftSet )\n\n");
		cat("Apologies for the change, but the new way is much better...\n")
		stop();
		
	}

	if (class(prior) != 'branchprior'){
		stop("object bprior is not of class branchprior\n");
	}
	
	dsc <- distinctShiftConfigurations(ephy, prior, BFcriterion);
	cfreq <- cumsum(dsc$frequency);
	cut <- min(which(cfreq >= set.limit));
	nodeset <- NULL;
 	
 	shiftnodes <- dsc$shifts[1:cut];
	indices <- dsc$samplesets[1:cut];
	frequency <- dsc$frequency[1:cut];
	cumulative <- cumsum(dsc$frequency)[1:cut];
	
	ephy$marg.probs <- dsc$marg.probs;
 	
 	ephy$shiftnodes <- shiftnodes;
 	ephy$indices <- indices;
 	ephy$frequency <- frequency;
 	ephy$cumulative <- cumulative;
 	ephy$coreshifts <- dsc$coreshifts;
 	ephy$BFcriterion <- BFcriterion;
 	ephy$set.limit <- set.limit;
 	ephy$number.distinct <- length(indices);
 	
	class(ephy) <- 'credibleshiftset';
	return(ephy);	
	
}






