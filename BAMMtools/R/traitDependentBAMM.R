
# Arguments:
#	ephy: a bammdata object
# rate: either "speciation", "extinction" or "net diversification". default to "speciation".
#	traits: a vector of trait data, with names corresponding to tips in the ephy object
#	reps: number of permutations to do
#	return.full: include the permutated and observed correlations in the retunred object?
#	method: either spearman,  pearson, mann-whitney, kruskal
#	logrates:  log-transform the rates before analysis? This can only matter for the pearson correlation
#	two.tailed: perform a two tailed statistical test? (in which case it will double the p-value)
# traitorder: for one tail test, specify the direction of correlation ("positive" or "negative") for countinues trait or a string 
#             indicating states of binary trait with increasing speciation rate, separted by comma (e.g., 'A, B'). Currently, this function only
#             perform two-tailed test for categorical data with more than two states.
# Returns:
#	estimate: the mean observed correlation coefficient
#	p.value: the probability that the observed correlation is less than or equal to
#				a sample from the null distribution. If you are doing a one-tailed
#				test on a continuous trait, a small pvalue means that your observed correlation is larger
#				than the null distribution. P-values approaching 1 mean that the observed
#				correlation is more negative than expected.
#	gen: if return.full=TRUE, the vector of generations sampled from the bammdata objected for permutation
#	obs.corr: if return.full=TRUE, a vector of observed correlation coefficients between traits and tip speciation rate for every sampled generation
#	null: if return.full=TRUE, a vector of permuted correlation coefficients for everty sample generation

# ephy<-ed
# traits<-bitrait
# reps<-100
# return.full=T
# method="m"
# logrates=T
# two.tailed=T
# nthreads=6
# traitorder='0,1';

##' @title STRAPP: STructured Rate Permutations on Phylogenies
##'
##' @description Given a \code{bammdata} object and a vector of (continuous)
##'     trait data, assess whether the correlation between the trait and bamm
##'     estimated speciation, extinction or net diversification rate is
##'     significant using permutation. A set of posterior samples is randomly
##'     drawn from the \code{bammdata} object. If the trait is continuous,
##'     this function calculates the correlation coefficients between the
##'     trait and tip rates (observed correlation), as well as that with
##'     permuted rates for each posterior sample. In a one-tailed test for
##'     positive correlations, the reported p-value is the proportion of the
##'     posterior samples in which the observed correlation is larger than the
##'     correlations calculated with permuted rates. In a two-tailed test, the
##'     reported p-value is the proportion of the posterior samples in which
##'     the null correlation is as extreme as the correlation observed. If the
##'     trait is binary, the U statistic of the Mann-Whitney test is
##'     calculated instead of correlation coefficients to assess whether there
##'     is a significant difference in rate between the two trait states. For
##'     categorical traits with more than two states, the Kruskal-Wallis rank
##'     sum statistic is used.
##'
##' @param ephy An object of class \code{bammdata}.
##' @param traits A vector of trait data, with names corresponding to tips in
##'     the \code{bammdata} object. It can be numeric or categorical. 
##' @param reps An integer specifying the number of permutations (i.e., the
##'     number of posterior samples to randomly draw with replacement from the
##'     \code{bammdata} object).
##' @param rate A character string specifying which estimated rate from the
##'     \code{bammdata} object to use for testing correlation, must be one of
##'     'speciation', 'extinction', 'net diversification' or 'trait'. Defaults to
##'     'speciation'. You can specify just the initial letter. Ignored for 
##' 	trait event data.
##' @param return.full A logical. If \code{TRUE}, the list of posterior
##'     samples, the observed correlation for each posterior sample, and the
##'     null distribution will be included in the returned object. Defaults to
##'     \code{FALSE}.
##' @param method A character string, must be one of 'spearman', 'pearson',
##'     'mann-whitney', or 'kruskal'. Defaults to 'spearman'. You can specify
##'     just the initial letter.
##' @param logrates A logical. If \code{TRUE} log-transform the rates before
##'     analysis. Defaults to \code{TRUE}. This can only matter for the
##'     pearson correlation. 
##' @param two.tailed A logical, used for continuous trait data. If
##'     \code{TRUE}, perform a two-tailed statistical test (i.e., if the null
##'     distribution is symmetric, it is equivalent to doubling the p-value).
##'     Defaults to \code{TRUE}.  
##' @param traitorder A character string specifying the direction of
##'     correlation for the alternative hypothesis. For continuous traits, it
##'     must be either "positive" or "negative"; only the initial letter is
##'     needed. For binary traits, it must be a string indicating states with
##'     increasing rate under the alternative hypothesis, separated by comma
##'     (e.g., 'A, B'). One-tailed test for categorical data with more than
##'     two states is not supported.
##' @param nthreads Number of threads to use for parallelization of the
##'     function. The R package \code{parallel} must be loaded for
##'     \code{nthreads > 1}.
##'
##' @details Tip rates --trait, speciation, extinction, or net diversification
##'     rates-- are permuted in a way such that pairwise covariances in rates
##'     between species are maintained. That is, tips with the same
##'     \code{tipStates} still have the same rate after permutation. Posterior
##'     samples are randomly selected with replacement from the
##'     \code{bammdata} object, so \code{reps} could be smaller or larger than
##'     the total number of samples in the object.
##'
##'     This function expects that the bamm-data object and the trait data
##'     have the same taxon set. It may be necessary to subset the trait data
##'     and/or run \code{\link{subtreeBAMM}} on the \code{bamm-data} object in
##'     order to meet this requirement.
##'
##' @return A list with the following components:
##'     \itemize{
##'         \item estimate: A numeric value for continous trait data: the
##'             average observed correlation between tip rates and the trait
##'             across the posterior samples. For categorical traits, it is
##'             a list showing the median species-specific rates for each
##'             trait state.
##'         \item p.value: A numeric value. The probability that the observed
##'             correlation is less than or equal to a sample from the null
##'             distribution.
##'         \item method: A character string, as input.
##'         \item rate: A character string, as input.
##'         \item two.tailed: A logical, as input.
##'         \item gen: An integer vector, recording which posterior samples
##'             were selected. Only present when \code{return.full} is
##'             \code{TRUE}.
##'         \item obs.corr: A numeric vector, the observed correlation
##'             coefficents for each posterior sample. Only present when
##'             \code{return.full} is \code{TRUE}. For binary traits, centered
##'             U statistics (U - n1* n2/2; where n1 and n2 are the number of
##'             species in each state of the binary trait) is reported.
##'         \item null: A numeric vector. The null distribution of
##'             correlation coefficients (or centered U statistics for binary
##'             traits) from permutation. Only present when \code{return.full}
##'             is \code{TRUE}.
##'     }
##'
##' @author Dan Rabosky, Huateng Huang
##'
##' @seealso \code{\link{subtreeBAMM}}
##'
##' @references \url{http://bamm-project.org/}
##'
##'     Rabosky, D. L. and Huang, H., 2015. A Robust Semi-Parametric Test for
##'     Detecting Trait-Dependent Diversification. Systematic Biology 65:
##'     181-193.
##'
##'     Rabosky, D. L. 2014. Automatic detection of key innovations, rate
##'     shifts, and diversity-dependence on phylogenetic trees. PLoS ONE
##'     9:e89543.
##'
##'     Rabosky, D. L., F. Santini, J. T. Eastman, S. A. Smith, B. L.
##'     Sidlauskas, J. Chang, and M. E. Alfaro. 2013. Rates of speciation and
##'     morphological evolution are correlated across the largest vertebrate
##'     radiation. Nature Communications DOI: 10.1038/ncomms2958.
##'
##' @examples
##' # using a small subset of the fish data set (300 taxa) in Rabosky et al.
##' # 2013. Nat. Com. paper
##' data(fishes, events.fishes)
##' xx <- getEventData(phy = fishes, eventdata = events.fishes, 
##'       nsamples = 500, type = "diversification")
##' # traits.fishes is the trait -- body size
##' data(traits.fishes)
##' x <- traitDependentBAMM(ephy = xx, traits = traits.fishes, reps = 1000,
##'                         return.full = TRUE, method = 's', logrates = TRUE,
##'                         two.tailed = TRUE)
##' @aliases strapp
##' @keywords nonparametric
##' @export
traitDependentBAMM <- function(ephy, traits, reps, rate = 'speciation', return.full = FALSE, method = 'spearman', logrates = TRUE, two.tailed = TRUE, traitorder = NA, nthreads = 1) {
   
	if (nthreads > 1) {
		if (!"package:parallel" %in% search()) {
			stop("Please load package 'parallel' for using the multi-thread option\n");
		}
	}
	
	if (ephy$type == 'trait') {
		rate <- 'trait'
	}
	
	ratetype.option <- c("speciation", "extinction", "net diversification", "trait");
	ratetype <- grep(paste0("^", rate), ratetype.option, value = TRUE, ignore.case = TRUE);
	  
  if (length(ratetype) == 0) {
    stop("Rate must be one of 'speciation', 'extinction', or 'net diversification', only the initial letter is needed.\n")
  }
  
	if (ratetype == "net diversification" & logrates == TRUE) {
		cat("WARNING: Net diversification might be negative and logged rates would then produce NaNs.\n");
	}

	if (ratetype == 'trait' & ephy$type == 'diversification') {
		stop('Rate must be either speciation or net diversification if ephy is a diversification analysis.\n');
	}

	if (ratetype %in% c('speciation', 'net diversification') & ephy$type == 'trait') {
		stop('Rate must be trait if ephy is a trait analysis.\n');
	}
	
	# check if species in ephy and traits match
	if (!identical(sort(names(traits)), sort(ephy$tip.label))) {
		stop('Species names in the bamm-data object and trait data are not identical. You may want to run subtreeBAMM and subset the trait data in order to reduce the two datasets to a common taxon set.')
	}

	method.option <- c("spearman",  "pearson", "mann-whitney", "kruskal");
	method <- grep(paste0("^", method), method.option, ignore.case = TRUE, value = TRUE);
	if (length(method) == 0) {
		stop("method must be one of 'spearman', 'pearson', 'mann-whitney', or 'kruskal', only the initial letter is needed");
	}

	# check if the trait is right class

	if (method == 'spearman' | method == "pearson") {
		if (! is.numeric(traits)){
			cat(paste0("selected ", method, ", but the trait is not numeric, converted the trait into a numeric vector\n"));
			traits <- as.numeric(traits);
		}
	} else if (method == "mann-whitney"| method == "kruskal") {
      
		if (length(unique(traits[! is.na(traits)])) == 1) {
			stop(paste("selected ", method, ", but the trait only has one level\n", sep = ''));
		}
		if (method == "mann-whitney") {
			if (length(unique(traits[! is.na(traits)])) > 2) {
				stop(paste("selected ", method, ", but the trait has more than two levels\n", sep = ''));
			}
		}
	}

	#check if the traitorder is specified
	trait.state <- NA;
	if (!two.tailed) {
		if (anyNA(traitorder)) {
			stop("selected one-tail test, but traitorder is not specified\n");
		}
		if ( method == "kruskal") {
			stop(" currently one-tail test is only available for continous or binary trait");
		}
		if (method == 'spearman' | method == "pearson") {
			direction.option <- c("positive", "negative");
			direction <- grep(paste0("^", traitorder), direction.option, ignore.case = TRUE, value = TRUE);
			if (length(direction) == 0) {
				stop(" for one-tail test with continous trait, traitorder must be either 'positive' or 'negative', only the initial letter is needed");
			} else {
				cat(paste0("select one-tailed ", method, " test\nAlternative hypothesis: the trait is ", direction, "ly correlated with speciation rate\n"));
			}
		} else {
			traitorder <- gsub(" ", "", traitorder);
			trait.state <- as.character(unlist(strsplit(x = traitorder, split = ",")));
			if (length(trait.state) != 2) {
				stop("please specify the traitorder for binary trait:\nTwo states separated by comma, and the state that is expected to have lower speciation rate first\n");
			} else {
				cat(paste0("selected one-tail ", method, " test\nAlternative hypothesis: species with trait ", trait.state[2], " has higher speciation rate than those with trait ", trait.state[1], "\n"));
			}
			for (i in trait.state) {
				if (sum(traits == i) == 0) {
					stop(paste("no species with state ", i," \n", sep = ''));
				}
			}
		}
	}
	if (ratetype == 'speciation') {
		tiprates <- ephy$tipLambda;
	} else if (ratetype == "extinction") {
		tiprates <- ephy$tipMu;
	} else if (ratetype == 'net diversification') {
		tiprates <- lapply(1:length(ephy$tipLambda), function(i) ephy$tipLambda[[i]] - ephy$tipMu[[i]]);
	} else if (ratetype == 'trait') {
		tiprates <- ephy$tipLambda;
	}
	tipstates <- ephy$tipStates;
	#tiprates <- tiprates[ephy$tip.label];
	traits <- traits[ephy$tip.label];
	stat.mu <- 0;
	if (method == "mann-whitney") {
		trait.stat.count <- table(traits);
		trait.stat.count <- trait.stat.count[! is.na(names(trait.stat.count))];
		stat.mu <- prod(trait.stat.count) / 2;
	}

	if (logrates) {
		tiprates <- lapply(1:length(tiprates), function(x) log(tiprates[[x]]));
	}

	#randomly sample generations from BAMM posterior
	gen <- sample(1:length(tiprates), size = reps, replace = TRUE);

	gen.tiprates <- list();
	for (l in 1:length(gen)) {
		gen.tiprates[[l]] <- data.frame(rates = tiprates[[gen[l]]], states = tipstates[[gen[l]]], stringsAsFactors = FALSE);
	}
  
	rm("tiprates", "tipstates");
	permute_tiprates <- function(m) {
		tt <- m$states;
		tlam <- m$rates;
		index <- unique(tt);
		lvec <- numeric(length(index));
		for (k in 1:length(index)) {
			lvec[k] <- tlam[tt == index[k]][1];
		}
		new_index <- sample(index, size = length(index), replace = FALSE);
		x <- rep(0,length(tt));
		for (xx in 1:length(index)) {
			x[which(tt == index[xx])] <- lvec[which(index == new_index[xx])];
		}
		x;   
	}
	if (nthreads > 1) {
		cl <- parallel::makePSOCKcluster(nthreads);
		p.gen.tiprates <- parallel::parLapply(cl, gen.tiprates, permute_tiprates);
		parallel::stopCluster(cl);
	} else {
		p.gen.tiprates <- lapply(gen.tiprates, permute_tiprates);
	}
	xgen.tiprates <- list();
	for (l in 1:length(gen)) {
		xgen.tiprates[[l]] <- gen.tiprates[[l]]$rates;
	}
	gen.tiprates <- xgen.tiprates; rm("xgen.tiprates");

	cortest <- function(rates, traits, method) {
	    if (sd(rates, na.rm = TRUE) == 0) {
	    	return(0);
	    } else {
			return(cor.test(rates, traits, method = method, exact = FALSE)$estimate);  
	    }
	}
	manntest <- function(rates, traits, two.tailed, trait.state) {
		if (two.tailed) {
			return(wilcox.test(rates ~ traits, exact = FALSE)$statistic);
		} else {
			return(wilcox.test(rates[which(traits == trait.state[2])], rates[which(traits == trait.state[1])], exact = FALSE)$statistic);
		}
	}

	kruskaltest <- function(rates, traits) {
		testres <- kruskal.test(rates ~ traits);
		# If there is no variation in rates (chi-squared value is NaN), then return a chi-squared value that has a P-value of 0.999, using the appropriate degrees of freedom
		if (is.na(testres$statistic)) {
			return(qchisq(p = 0.999, df = testres$parameter));
		} else {
			return(testres$statistic);
		}
	}
 

	if (nthreads > 1) {
		cl <- parallel::makePSOCKcluster(nthreads);
		if (method == 'spearman' | method == "pearson") {
			obs <- parallel::parLapply(cl, gen.tiprates,cortest, traits, method);
			permu <- parallel::parLapply(cl, p.gen.tiprates, cortest, traits, method);
		} else if (method == "mann-whitney") {
			obs <- parallel::parLapply(cl, gen.tiprates, manntest, traits, two.tailed, trait.state);
			permu <- parallel::parLapply(cl, p.gen.tiprates, manntest, traits, two.tailed, trait.state);
		} else {
			obs <- parallel::parLapply(cl, gen.tiprates, kruskaltest, traits);
			permu <- parallel::parLapply(cl, p.gen.tiprates, kruskaltest, traits);
		}
		parallel::stopCluster(cl);
	} else {
		if (method == 'spearman' | method == "pearson") {
			obs <- lapply(gen.tiprates, cortest, traits, method);
			permu <- lapply(p.gen.tiprates, cortest, traits, method);
		} else if (method == "mann-whitney") {      
			obs <- lapply(gen.tiprates, manntest, traits, two.tailed, trait.state);
			permu <- lapply(p.gen.tiprates, manntest, traits, two.tailed, trait.state);
		} else {    
			obs <- lapply(gen.tiprates, kruskaltest, traits);
			permu <- lapply(p.gen.tiprates, kruskaltest, traits);    
		}
	}
	obs <- unlist(obs);
	permu <- unlist(permu);
	obs <- obs - stat.mu;
	permu <- permu - stat.mu;

  

	if (two.tailed) {
		pval <- sum(abs(obs) <= abs(permu)) / length(permu);
	} else {
		if (method == "spearman" | method == "pearson") {
			if (direction == 'positive') {
				pval <- sum(obs <= permu) / length(permu);
			} else {
				pval <- sum(obs >= permu) / length(permu);
			}
		} else {
			pval <- sum(obs <= permu) / length(permu);
		}
	}
	if (method == "spearman" | method == "pearson") {
		obj <- list(estimate = mean(as.numeric(obs)), p.value = pval, method = method, two.tailed = two.tailed);
	} else {
		if (ratetype == 'speciation') {
			ave.tiprate <- getTipRates(ephy)$lambda.avg;
		} else if (ratetype == 'extinction') {
			ave.tiprate <- getTipRates(ephy)$mu.avg;
		} else if (ratetype == 'net diversification') {
			ave.tiprate <- getTipRates(ephy)$lambda.avg - getTipRates(ephy)$mu.avg;
		} else if (ratetype == 'trait') {
			ave.tiprate <- getTipRates(ephy)$beta.avg;
		}
		l <- lapply(unique(traits[!is.na(traits)]), function(x) {
			median(ave.tiprate[which(traits == x)], na.rm = TRUE);
		});
		names(l) <- as.character(unique(traits[! is.na(traits)]));
		obj <- list(estimate = l, p.value = pval, method = method, two.tailed = two.tailed);
	}
	obj$rate <- ratetype;
	if (return.full) {
		obj$obs.corr <- as.numeric(obs);
		obj$gen <- gen;
		obj$null <- as.numeric(permu);
	}

	return(obj);
}

