
# phy is a phylogenetic tree
# traits is a file name to your BAMM formatted data.
# total.taxa is only if you have incomplete sampling
# 
# 

##' @title Set BAMM Priors
##'
##' @description Set priors for \code{BAMM} analysis.
##'
##' @param phy An object of class \code{phylo}, e.g., the phylogenetic tree
##'     you will analyze with \code{BAMM}.
##' @param total.taxa If doing speciation-extinction analysis, the total
##'     number of taxa in your clade. If your tree contains all taxa in the
##'     clade (100\% sampling), then leave this as \code{NULL}.
##' @param traits A filename where the trait data (\code{BAMM} format) are
##'     stored, or a numeric vector named according to the tips in \code{phy}.
##'     Leave \code{NULL} if doing a speciation-extinction analysis.
##' @param outfile Filename for outputting the sample priors block. If
##'     \code{NULL}, then a vector is returned instead.
##' @param Nmax If analyzing a very large tree for phenotypic evolution, uses
##'     only this many taxa to estimate priors for your dataset. Avoid matrix
##'     inversion issues with large numbers of tips. 
##' @param suppressWarning Logical. If \code{TRUE}, then the warning about
##'     setting the Poisson rate prior is suppressed. Only applies if
##'     \code{outfile = NULL}.
##'
##' @details This is a "quick and dirty" tool for identifying approximately
##'     acceptable priors for a \code{BAMM} analysis. We have found that
##'     choice of prior can have a substantial impact on \code{BAMM} analyses.
##'     It is difficult to simply set a default prior that applies across
##'     datasets, because users often have trees with branch lengths in very
##'     different units (e.g., numbers of substitutions versus millions of
##'     years). Hence, without some careful attention, you can inadvertently
##'     specify some very bad prior distributions. This function is designed
##'     to at least put you in the right ballpark for decent prior
##'     distributions, but there are no guarantees that these are most
##'     appropriate for your data.
##'
##'     The general rules applied here are: 
##'
##'     For the \code{lambdaInitPrior}, we estimate the speciation rate of the
##'     data under a pure birth model. We then set this prior to give an
##'     exponential distribution with a mean five times greater than this
##'     computed pure birth speciation estimate.
##'
##'     The \code{lambdaShiftPrior} is the standard deviation of the normal
##'     prior on the exponential change parameter k. We set the prior
##'     distribution based on the age of the root of the tree. We set the
##'     standard deviation of this distribution such that 2 standard
##'     deviations give a parameter that will yield a 90\% decline in the
##'     initial speciation rate between the root of the tree and the tips.
##'     The basic model is lambda(t) = lambda_0 * exp(k * t). This is a
##'     straightforward calculation: let x = -log(0.1) / TMAX, where TMAX is
##'     the age of the tree. Then set the standard deviation equal to (x / 2).
##'
##'     We set \code{muInitPrior} equal to \code{lambdaInitPrior}.  
##'
##'     For trait evolution, we first compute the maximum likelihood estimate
##'     of the variance of a Brownian motion process of trait evolution. The
##'     prior \code{betaInitPrior} is then set to an exponential distribution
##'     with a mean 5 times greater than this value (similar to what is done
##'     for lambda and mu, above).  
##'
##'     This function generates an output file containing a prior parameters
##'     block that can be pasted directly into the priors section of your
##'     \code{BAMM} control file.
##'
##' @return The function does not return anything. It simply performs some
##'     calculations and writes formatted output to a file. However, if
##'     \code{outfile = NULL}, then a named vector is returned.
##'
##' @author Dan Rabosky
##'
##' @examples
##' # for diversification analyses
##' data(whales)
##' setBAMMpriors(phy = whales, total.taxa = 89, outfile = NULL)
##'
##' # for trait analyses
##' data("primates")
##' data("mass.primates")
##'
##' ## create a named vector of trait values
##' mass <- setNames(mass.primates[,2], mass.primates[,1])
##'
##' setBAMMpriors(phy = primates, traits = mass, outfile = NULL)
##' @keywords models
##' @export
setBAMMpriors <- function(phy, total.taxa = NULL, traits=NULL, outfile = 'myPriors.txt', Nmax = 1000, suppressWarning = FALSE){
	
	if (is.ultrametric(phy)) {
		mbt <- max(branching.times(phy));
	} else {
		mbt <- max(NU.branching.times(phy));
	}
	if (is.null(total.taxa)){
		total.taxa <- length(phy$tip.label);
	}
	if (!is.null(total.taxa) & total.taxa <= 1) {
		stop("total.taxa is a count, not a percent.");
	}
	
	if (length(phy$tip.label) > total.taxa) {
		stop("Your tree has more tips than total.taxa...");
	}
	
	if (is.null(traits)){
		
		pb <- (log(total.taxa) - log(2)) / mbt;
		lamprior <- 1 / (pb * 5);
		lamrootprior <- 1 / (pb * 1);
		k1 <- log(0.1) / mbt;
		kprior <- -1 * (k1 / 2);	
		
		s1 <- '###############################################';
		s2 <- '# Prior block chosen by BAMMtools::setBAMMpriors';
		s3 <- 'expectedNumberOfShifts = 1.0';
		s4 <- paste('lambdaInitPrior = ', lamprior, sep='');
		s5 <- paste('lambdaShiftPrior = ', kprior, sep='');
		s6 <- paste('muInitPrior = ', lamprior, sep='');
		s7 <- paste('#### End Prior block\n######################\n\n');
		ss <- paste(s1,s2,s3,s4,s5,s6,s7, sep='\n\n');
		if (!is.null(outfile)) {
			write(ss, file = outfile, sep='');
		} else {
			res <- as.data.frame(cbind(c('expectedNumberOfShifts', 'lambdaInitPrior', 'lambdaShiftPrior', 'muInitPrior'), c(1.0, lamprior, kprior, lamprior)), stringsAsFactors=FALSE);
			res[,2] <- as.numeric(res[,2]);
			colnames(res) <- c('param','value');
		}
	} else {
		if (is.character(traits)) {
			x <- read.table(file = traits, sep = '\t', stringsAsFactors = FALSE, header = FALSE);
			tvec <- x[,2];
			names(tvec) <- x[,1];
		} else {
			tvec <- traits;
		}
		not.in.tree <- setdiff(phy$tip.label, names(tvec));
		not.in.traits <- setdiff(names(tvec), phy$tip.label);
		bad <- c(not.in.tree, not.in.traits);
		if (length(bad) > 0){
			cat('Taxon names do not match between tree and trait dataset\n');
			cat('Printing mismatched names:\n\n');
			for (i in bad){
				cat(i, '\n');
			}
			stop('Names in trait dataset must match those in tree\n');
		}
		bad <- which(is.na(tvec)) #check for missing data
		if (length(bad) > 0) {
			cat('\nSome taxa are missing data:\n\n');
			for (i in bad) {
				cat(names(tvec)[i], '\n');
			}
			stop('All taxa must have trait data.');
		}
		tvec <- tvec[phy$tip.label];
		if (length(phy$tip.label) > Nmax){
			ss <- sample(phy$tip.label, size=Nmax);
			drop <- setdiff(phy$tip.label, ss);
			phy <- drop.tip(phy, tip = drop);
			tvec <- tvec[phy$tip.label];
		}
		pmean <- phylogeneticMean(tvec, phy)$beta;
		
		betaprior <- 1/(pmean * 5);
		betarootprior <- 1/(pmean * 1);		
		
		k1 <- log(0.1) / mbt;
		kprior <- -1 * (k1 / 2);	
		
		s1 <- '###############################################';
		s2 <- '# Prior block chosen by BAMMtools::setBAMMpriors';
		s3 <- 'expectedNumberOfShifts = 1.0';
		s4 <- paste('betaInitPrior = ', betaprior, sep='');
		s5 <- paste('betaShiftPrior = ', kprior, sep='');
		s6 <- paste('useObservedMinMaxAsTraitPriors = 1');
		s7 <- paste('#### End Prior block\n######################\n\n');
		ss <- paste(s1,s2,s3,s4,s5,s6,s7,sep='\n\n');
		if (!is.null(outfile)) {
			write(ss, file = outfile, sep='');
		} else {
			res <- as.data.frame(cbind(c('expectedNumberOfShifts', 'betaInitPrior', 'betaShiftPrior', 'useObservedMinMaxAsTraitPriors'), c(1.0, betaprior, kprior, 1)), stringsAsFactors=FALSE);
			res[,2] <- as.numeric(res[,2]);
			colnames(res) <- c('param','value');
		}
				
	}

	if (!is.null(outfile)) {	
		cat('\nPrior block written to file << ', outfile, " >>\n", sep='');
		cat('Copy and paste the contents of the file into the\n');
		cat('priors block of your BAMM input file\n');
	}

	if (!suppressWarning & !is.null(outfile)) {
		cat('\nThis function simply sets the expectedNumberOfShifts to 1;\n');
		cat('This is a parameter you may need to vary to achieve good convergence\n');
		cat('with your data.\n');
	}

	if (is.null(outfile)) {
		res <- setNames(res[,2], res[,1])
		return(res);
	}	
}







