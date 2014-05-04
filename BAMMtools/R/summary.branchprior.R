summary.branchprior <- function(object, ...){
	if (class(object) != 'branchprior'){
		stop('obj x is not of class branchprior');
	}
	per <- round(object$criterion*100, 2);
	invper <- 100 - per;
	cat('\nObject of class "branchprior".\n');
	cat("Branch lengths are equal to the prior probability of a\n", sep='');
	cat("shift event occurring on the focal branch\n");
}

