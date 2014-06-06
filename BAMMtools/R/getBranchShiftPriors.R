
getBranchShiftPriors <- function(phy, priordata){
	
	if (class(priordata) == 'character'){
		prior <- read.csv(priordata, header=T);
	}else if (class(priordata) == 'data.frame'){
		prior <- priordata;

	}else{
		stop("invalid priordata argument (wrong class) in getBranchShiftPriors\n");
	}
	
	if (sum(c('N_shifts', 'prob') %in% colnames(prior)) != 2){
		cat('Invalid data.frame passed to getBranchShiftPriors\n');
		cat('File must have columns with names << N_shifts>> and << prob >>\n');
		stop("Please check!");
	}
	
	
	prior <- prior[prior$N_shifts > 0 & prior$prob > 0, ];
	
	pvec <- phy$edge.length / sum(phy$edge.length);
	
	pp <- numeric(length(phy$edge.length));
	
	
	
	for (i in 1:nrow(prior)){
		# probability of getting 0 shifts on branch given ns total 
		#  given the branch lengths etc
		#	weighted by the probability of that shift category
		
		pp <- pp + (1 - dbinom(0, prior$N_shifts[i], prob=pvec))*prior$prob[i];
		
	}
	
	obj <- phy;	
	obj$edge.length <- pp;
	class(obj) <- 'branchprior';
	return(obj);
}

