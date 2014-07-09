
computeBayesFactors <- function(postdata, priordata, burnin = 0.1, modelset = NULL, ...){

	if (hasArg("strict") | hasArg("threshpost") | hasArg("threshprior") | hasArg("nbprior")){
 		cat("Error - you have specified some argument names that have been deprecated\n");
 		cat("in this version of BAMMtools. Check the help file on this function\n");
 		cat("to see what has changed\n\n");
		stop();
		
	}


	if (class(postdata) == 'character'){
		dpost <- read.csv(postdata, header=T);
	}else if (class(postdata) == 'data.frame'){
		dpost <- postdata;
	}else{
		stop("invalid postdata argument (wrong class) in computeBayesFactors\n");
	}

	if (class(priordata) == 'character'){
		prior <- read.csv(priordata, header=T);
	}else if (class(priordata) == 'data.frame'){
		prior <- priordata;
	}else{
		stop("invalid priordata argument (wrong class) in computeBayesFactors\n");
	}
 
	dpost <- dpost[floor(burnin*nrow(dpost)):nrow(dpost), ];
 
	tx <- table(dpost$N_shifts) / nrow(dpost);
	
	post <- data.frame(N_shifts=as.numeric(names(tx)), prob=as.numeric(tx));
	
 	### Restrict to set of models that have been sampled:
 	
 	inboth <- intersect(prior$N_shifts, post$N_shifts);
 
 	if (length(inboth) == 0){
 		cat("No overlap between models sampled during simulation of the prior\n");
 		cat("and those sampled during simulation of posterior\n");
		cat("Bayes factors cannot be computed in this case due to difficulty\n");
		cat("approximating model probabilities for (very) rarely sampled models\n");
		stop();
 	}else{
 		prior <- prior[prior$N_shifts %in% inboth, ];
 		post <- post[post$N_shifts %in% inboth, ];
 	}
	
	not.in.obs <- setdiff(modelset, inboth);
	
	if (!is.null(modelset)){
		inboth <- intersect(inboth, modelset);
		if (length(inboth) == 0){
			cat("No overlap between sampled models and those for which\n");
			cat("you wish to compute Bayes factors\n");
			stop();
		}else if (length(not.in.obs) > 0){
			cat("\nThe modelset you specified includes models that were not sampled by BAMM\n");
			cat("Consequently, you cannot estimate their probabilities with any degree\n");
			cat("of accuracy. They will be excluded from the Bayes factor matrix.\n\n");
		}
	} 
 
	mm <- matrix(NA, nrow=length(inboth), ncol=length(inboth));
	rownames(mm) <- inboth;
	colnames(mm) <- inboth;
 
	for (i in 1:length(inboth)){
			mi <- inboth[i];
		for (j in 1:length(inboth)){
			
			mj <- inboth[j];
				
			prior_odds <- prior$prob[prior$N_shifts == mi] / prior$prob[prior$N_shifts == mj];
			post_odds <- post$prob[post$N_shifts == mi] / post$prob[post$N_shifts == mj];
			
			mm[i,j] <- post_odds * (1 / prior_odds);
			
		}	
		
	}
	
	return(mm);
	
}



