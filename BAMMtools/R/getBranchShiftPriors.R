
getBranchShiftPriors <- function(phy, priordata){
	
	if (class(priordata) == 'character'){
		prior <- read.csv(priordata, header=T);
	}else if (class(priordata) == 'data.frame'){
		prior <- priordata;
	}else{
		stop("invalid priordata argument (wrong class) in getBranchShiftPriors\n");
	}
	
	tx <- table(prior$N_shifts) / nrow(prior);
	
	tx <- tx[names(tx) != '0'];
	ns <- as.numeric(names(tx));
	
	pvec <- phy$edge.length / sum(phy$edge.length);
	
	pp <- numeric(length(phy$edge.length));
	
	for (i in 1:length(tx)){
		
		pp <- pp + (1 - dbinom(0, ns[i], prob=pvec))*tx[i];
		
	}
	
	obj <- phy;	
	obj$edge.length <- pp;
	class(obj) <- 'branchprior';
	return(obj);
}

