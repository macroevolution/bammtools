##' @title Create \code{bammdata} object from MCMC output
##'
##' @description \code{getEventData} Reads shift configuration data (the
##'     "event data" output) from a \code{BAMM} analysis and creates a
##'     \code{bammdata} object. The \code{bammdata} object is fundamental
##'     for extracting information about macroevolutionary rate variation
##'     through time and among lineages.
##'
##' @param phy An object of class \code{phylo} - specifically, the
##'     time-calibrated tree that was analyzed with \code{BAMM}.
##'     Alternatively, a character string specifying the path to a
##'     newick-formatted tree.
##' @param eventdata A character string specifying the path to a \code{BAMM}
##'     event-data file. Alternatively, an object of class \code{data.frame}
##'     that includes the event data from a \code{BAMM} run.
##' @param burnin A numeric indicating the fraction of posterior samples to
##'     discard as burn-in.
##' @param nsamples An integer indicating the number of posterior samples to
##'     include in the \code{bammdata} object. May be \code{NULL}.
##' @param verbose A logical. If \code{TRUE} progess is outputted to the
##'     console. Defaults to \code{FALSE}.
##' @param type A character string. Either "diversification" or "trait"
##'     depending on your \code{BAMM} analysis.
##'
##' @details In the \code{BAMM} framework, an "event" defines a
##'     macroevolutionary process of diversification or trait evolution. Every
##'     sample from the posterior includes at least one process, defined by
##'     such an "event". If a given sample includes just a single event, then
##'     the dynamics of diversification or trait evolution can be described
##'     entirely by a single time-constant or time-varying process that begins
##'     at the root of the tree. Any sample from the posterior distribution
##'     may include a complex mixture of distinct processes. To represent
##'     temporal heterogeneity in macroevolutionary rates, \code{BAMM} models
##'     a rate \eqn{R}, e.g. speciation, as a function that changes
##'     exponentially with time:
##'
##'     \eqn{R(t) = R(0)*exp(b*t)}.
##'
##'     Here \eqn{R(0)} is the initial rate and \eqn{b} is a parameter
##'     determining how quickly that rate grows or decays with time. 
##'
##'     The \code{eventdata} file (or data frame) is a record of events and
##'     associated parameters that were sampled with \code{BAMM} during
##'     simulation of the posterior with reversible jump MCMC. This complex,
##'     information-rich file is processed into a \code{bammdata} object,
##'     which serves as the core data object for numerous downstream analyses.
##'     From a \code{bammdata} object, you can summarize rate variation
##'     through time, among clades, extract locations of rate shifts,
##'     summarize clade-specific rates of speciation and extinction, and more.
##'
##'     In general, the user does not need to be concerned with the details of
##'     a \code{bammdata} object. The object is used as input by a number of
##'     \code{BAMMtools} functions. 
##'
##'     The parameter \code{nsamples} can be used to reduce the total amount
##'     of data included in the raw eventdata output from a \code{BAMM} run.
##'     The final \code{bammdata} object will consist of all data for
##'     \code{nsamples} from the posterior. These \code{nsamples} are equally
##'     spaced after discarding some \code{burnin} fraction as "burn-in". If
##'     \code{nsamples} is set to \code{NULL}, the \code{bammdata} object will
##'     include all samples in the posterior after discarding the
##'     \code{burnin} fraction.
##'
##' @return A list with many components:
##' \itemize{
##'     \item{edge} {See documentation for class \code{phylo} in package ape.}
##'     \item{Nnode} {See documentation for class \code{phylo} in package
##'         ape.}
##'     \item{tip.label} {See documentation for class \code{phylo} in package
##'         ape.}
##'     \item{edge.length} {See documentation for class \code{phylo} in
##'         package ape.}
##'     \item{begin} {The beginning time of each branch in absolute time (the
##'         root is set to time zero)}
##'     \item{end} {The ending time of each branch in absolute time.}
##'     \item{numberEvents} {An integer vector with the number of events
##'         contained in \code{phy} for each posterior sample. The length of
##'         this vector is equal to the number of posterior samples in the
##'         \code{bammdata} object.}
##'     \item{eventData} {A list of dataframes. Each element is a single
##'         posterior sample. Each row in a dataframe holds the data for a
##'         single event. Data associated with an event are: \code{node} - a
##'         node number. This identifies the branch where the event
##'         originates. \code{time} - this is the absolute time on that branch
##'         where the event originates (with the root at time 0). \code{lam1}
##'         - an initial rate of speciation or trait evolution. \code{lam2} -
##'         a decay/growth parameter. \code{mu1} - an initial rate of
##'         extinction. \code{mu2} - a decay/growth parameter. \code{index} -
##'         a unique integer associated with the event. See 'Details'.}
##'     \item{eventVectors} {A list of integer vectors. Each element is a
##'         single posterior sample. For each branch in \code{phy} the index
##'         of the event that occurs along that branch. Branches are ordered
##'         increasing here and elsewhere.}
##'     \item{eventBranchSegs} {A list of matrices. Each element is a single
##'         posterior sample. Each matrix has four columns: \code{Column 1}
##'         identifies a node in \code{phy}. \code{Column 2} identifies the
##'         beginning time of the branch or segment of the branch that
##'         subtends the node in \code{Column 1}. \code{Column 3} identifies
##'         the ending time of the branch or segment of the branch that
##'         subtends the node in \code{Column 1}. \code{Column 4} identifies
##'         the index of the event that occurs along the branch or segment of
##'         the branch that subtends the node in \code{Column 1}.}
##'     \item{tipStates} {A list of integer vectors. Each element is a single
##'         posterior sample. For each tip the index of the event that occurs
##'         along the branch subtending the tip. Tips are ordered increasing
##'         here and elsewhere.}
##'     \item{tipLambda} {A list of numeric vectors. Each element is a single
##'         posterior sample. For each tip the rate of speciation or trait
##'         evolution at the end of the terminal branch subtending that tip.}
##'     \item{tipMu} {A list of numeric vectors. Each element is a single
##'         posterior sample. For each tip the rate of extinction at the end
##'         of the terminal branch subtending that tip. Meaningless if working
##'         with \code{BAMM} trait results.}
##'     \item{meanTipLambda} {For each tip the mean of the marginal posterior
##'         density of the rate of speciation or trait evolution at the end of
##'         the terminal branch subtending that tip.}
##'     \item{meanTipMu} {For each tip the mean of the marginal posterior
##'         density of the rate of extinction at the end of the terminal
##'         branch subtending that tip. Meaningless if working with
##'         \code{BAMM} trait results.}
##'     \item{type} {A character string. Either "diversification" or "trait"
##'         depending on your \code{BAMM} analysis.}
##'     \item{downseq} {An integer vector holding the nodes of \code{phy}. The
##'         order corresponds to the order in which nodes are visited by a
##'         pre-order tree traversal.}
##'     \item{lastvisit} {An integer vector giving the index of the last node
##'         visited by the node in the corresponding position in
##'         \code{downseq}. \code{downseq} and \code{lastvisit} can be used to
##'         quickly retrieve the descendants of any node. e.g. the descendants
##'         of node 89 can be found by
##'         \code{downseq[which(downseq==89):which(downseq==lastvisit[89])}.}
##' }
##'
##' @note Currently the function does not check for duplicate tip labels in
##'     \code{phy}, which may cause the function to choke.
##'
##' @author Dan Rabosky, Mike Grundler
##'
##' @seealso \code{\link{summary.bammdata}}, \code{\link{plot.bammdata}},
##'     \code{\link{dtRates}}.
##' 
##' @references \url{http://bamm-project.org/}
##'
##' @examples
##' data(primates, events.primates)
##' xx <- getEventData(primates, events.primates, burnin=0.25, nsamples=500,
##'                    type = 'trait')
##' 
##' # compute mean phenotypic rate for primate body size evolution:
##' brates <- getCladeRates(xx)
##' mean(brates$beta)
##' 
##' # Plot rates:
##' plot(xx)
##' @keywords models
##' @export
getEventData <- function(phy, eventdata, burnin=0, nsamples = NULL, verbose=FALSE, type = 'diversification')
{	
	if (type != 'diversification' & type != 'trait') {
		stop("Invalid 'type' specification. Should be 'diversification' or 'trait'");
	}
	
	if (inherits(phy, 'character')) {
		phy <- read.tree(phy);
	}
	
	phy$node.label <- NULL;
	
	if (any(is.null(c(phy$begin, phy$end)))) {
		phy <- getStartStopTimes(phy)
	}
	
# Getting branching times direct from
#		the begin and end components of phy
#		should be able to now handle non-ultrametric trees.
	
	maxbt <- max(phy$end)
	nodes <- (length(phy$tip.label) + 1):(2*length(phy$tip.label) - 1)
	bt <- numeric(length(nodes))
	names(bt) <- nodes
	for (i in 1:length(bt)){
		tt <- phy$begin[phy$edge[,1] == nodes[i]][1]
		bt[i] <- maxbt - tt
	}
	
	
	########
	
	if (inherits(eventdata, 'data.frame')) {
		cat("Processing event data from data.frame\n");
		uniquegens <- sort(unique(eventdata[,1]));
	} 
	else if (inherits(eventdata, 'character')) {
		cat("Reading event datafile: ", eventdata, "\n\t\t...........");
		eventdata <- read.csv(eventdata, header=TRUE, stringsAsFactors=FALSE);
 		uniquegens <- sort(unique(eventdata[,1]));
 		cat("\nRead a total of", length(uniquegens), "samples from posterior\n");				
	} 
	else {
		err.string = c('eventdata arg invalid\n\nType is ', class(eventdata), '\n', sep='');
		stop(err.string);
	}
 	
 	samplestart <- uniquegens[floor(burnin*length(uniquegens))];
 	if (!length(samplestart)) {
 		samplestart <- 0;
 	}
 	uniquegens <- uniquegens[uniquegens >= samplestart];
 
 	if (is.null(nsamples)) {
 		nsamples <- length(uniquegens);
 	} 
 	else if (nsamples > length(uniquegens)) {
 		nsamples <- length(uniquegens);
 	}
	
	eventVectors <- vector("list",nsamples);
	eventData <- vector("list",nsamples);
	tipStates <- vector("list",nsamples);
	eventBranchSegs <- vector("list",nsamples);
	
	tipLambda <- vector("list",nsamples);
	tipMu <- vector("list",nsamples);
	
	goodsamples <- uniquegens[seq.int(1, length(uniquegens), length.out=nsamples)];
 	 
 	cat('\nDiscarded as burnin: GENERATIONS < ', goodsamples[1]);
 	cat("\nAnalyzing ", length(goodsamples), " samples from posterior\n");

 	numberEvents <- length(goodsamples); # vector to hold number of events per sample
 
 	cat('\nSetting recursive sequence on tree...\n');
 	phy <- getRecursiveSequence(phy);
 	cat('\nDone with recursive sequence\n\n');
 
	######### Get ancestors for each pair of taxa
	if (verbose) {
		cat("Start preprocessing MRCA pairs....\n");
	}	
		
	x2 <- eventdata[eventdata$generation %in% goodsamples, ];
		
	uniquePairSet <- matrix(NA, nrow=nrow(x2), ncol=2);	
	uniquePairNode <- numeric(nrow(x2));
	
	uniquePairSet[,1] <- as.integer(match(x2$leftchild, phy$tip.label));
	uniquePairSet[,2] <- as.integer(match(x2$rightchild, phy$tip.label, nomatch = 0L));
	uniquePairNode <- getmrca(phy, uniquePairSet[,1],uniquePairSet[,2]);
	
	if (verbose) {
		cat("Done preprocessing MRCA pairs....\n");
	}	

	####### Done with risky sstuff
	
	meanTipMu <- numeric(length(phy$tip.label));
	
 	meanTipLambda <- numeric(length(phy$tip.label)); 
 
 	stoptime <- maxbt;
 	
 	for (i in 1:length(goodsamples)) {
  		tmpEvents <- x2[x2[,1] == goodsamples[i], ];
		
		if (verbose) cat('Processing event: ', i, '\n');		
 
 		tm <- as.numeric(tmpEvents[,4]); # abs time of event
 		lam1 <- as.numeric(tmpEvents[,5]); # lambda parameter 1
 		lam2 <- as.numeric(tmpEvents[,6]); # lambda parameter 2
 		if (type == 'diversification') {	
			mu1 <- try(as.numeric(tmpEvents[, 7]),silent=TRUE); # mu parameter 1
			if (inherits(mu1,"try-error")) {
				stop("Unidentified column in event data file. Maybe try setting argument 'type = trait'");
			}
 			mu2 <- as.numeric(tmpEvents[, 8]); #mu parameter 2 
		} 
		else { #for bamm trait data we set the mu columns to zero because those params don't exist	
			mu1 <- rep(0, nrow(tmpEvents)); 
 			mu2 <- rep(0, nrow(tmpEvents)); 
		}	
		
 		# Get subtending node for each event:
 		nodeVec <- uniquePairNode[x2[,1] == goodsamples[i]];
 		
 		if (sum(nodeVec == 0)  > 0) {
			stop('Failed to assign event to node\n');
		}
		
		# make a dataframe:
		dftemp <- data.frame(node=nodeVec, time=tm, lam1=lam1, lam2=lam2, mu1=mu1, mu2=mu2, stringsAsFactors=FALSE);
		
		dftemp <- dftemp[order(dftemp$time), ];
		dftemp$index <- 1:nrow(dftemp);
		rownames(dftemp) <- NULL;
		
		statevec <- rep(1, nrow(phy$edge));

		if (nrow(dftemp) > 1) {
			for (k in 2:nrow(dftemp)) {
				s1 <- which(phy$downseq == dftemp$node[k]);
				s2 <- which(phy$downseq == phy$lastvisit[dftemp$node[k]]);
				descSet <- phy$downseq[s1:s2];
				isDescendantNode <- phy$edge[,2] %in% descSet;				
				statevec[isDescendantNode] <- k;
			}				
		}

 		tmpEventSegMat <- matrix(0, nrow=(max(phy$edge) + nrow(dftemp) - 2), ncol=4);
 		
		#non.root <- c(1:length(phy$tip.label), (length(phy$tip.label)+2):max(phy$edge));
		non.root <- c(1:length(phy$tip.label), unique(phy$edge[,1]));
		non.root <- non.root[-match(length(phy$tip.label)+1, non.root)];
		pos <- 1;	
		
		is_noEventBranch = !(phy$edge[,2] %in% dftemp$node);
		
		if (sum(is_noEventBranch) > 0) {
			
			tmpEventSegMat[1:sum(is_noEventBranch), 1] <- phy$edge[,2][is_noEventBranch];
			tmpEventSegMat[1:sum(is_noEventBranch), 2] <- phy$begin[is_noEventBranch];
	 		tmpEventSegMat[1:sum(is_noEventBranch), 3] <- phy$end[is_noEventBranch];
	 		tmpEventSegMat[1:sum(is_noEventBranch), 4] <- statevec[is_noEventBranch];
	 			
		} else {
			tempEventSegMat <- cbind(phy$edge[, 2], phy$begin, phy$end, statevec);
		}
 		
		eventnodeset <- intersect(non.root, dftemp$node);
		pos <- 1 + sum(is_noEventBranch);
		for (k in eventnodeset) {
			events.on.branch <- dftemp[dftemp$node == k, ];
			events.on.branch <- events.on.branch[order(events.on.branch$time), ];
				
			fBranch <- phy$edge[,2] == k;
 			start.times <- c(phy$begin[fBranch], events.on.branch$time);
			stop.times <- c(events.on.branch$time, phy$end[fBranch]);
			parent <- phy$edge[,1][phy$edge[,2] == k];
			if (parent == (length(phy$tip.label) + 1)) {
				# Parent is root:
				proc.set <- c(1, events.on.branch$index);	
			} else {
				proc.set <- c(statevec[phy$edge[,2] == parent], events.on.branch$index);			
			}
				
 			zzindex = pos:(pos + nrow(events.on.branch));	
				
			tmpEventSegMat[zzindex, 1] <- rep(k, length(zzindex));
			tmpEventSegMat[zzindex, 2] <- start.times;
			tmpEventSegMat[zzindex, 3] <- stop.times;
			tmpEventSegMat[zzindex, 4] <- proc.set;		
			pos <- pos + 1 + nrow(events.on.branch);
		}
		
 		tmpEventSegMat <- tmpEventSegMat[order(tmpEventSegMat[,1]), ];
 	
 		eventBranchSegs[[i]] <- tmpEventSegMat;

		tipstates <- numeric(length(phy$tip.label));
		tipstates <- statevec[phy$edge[,2] <= phy$Nnode + 1];
		tipstates <- tipstates[order(phy$edge[phy$edge[,2] <= phy$Nnode + 1, 2])];
		
 		### Compute tip rates:
		
		#tiplam <- dftemp$lam1[tipstates] * exp(dftemp$lam2[tipstates] * (stoptime - dftemp$time[tipstates]));
		tiplam <- exponentialRate(stoptime - dftemp$time[tipstates], dftemp$lam1[tipstates], dftemp$lam2[tipstates]);
		tipmu <- dftemp$mu1[tipstates];
		
		meanTipMu <- meanTipMu + tipmu/nsamples;
		meanTipLambda <- meanTipLambda + tiplam/nsamples;
		
		### List assignments and metadata across all events:
		eventData[[i]] <- dftemp;	
		eventVectors[[i]] <- statevec;
		numberEvents[i] <- nrow(dftemp);
		tipStates[[i]] <- tipstates;
		
		tipLambda[[i]] <- tiplam;
		tipMu[[i]] <- tipmu;	
 	}
 	
	phy$numberEvents <- numberEvents;
	phy$eventData <- eventData;
	phy$eventVectors <- eventVectors;
	phy$tipStates <- tipStates;
	phy$tipLambda <- tipLambda;
	phy$meanTipLambda <- meanTipLambda;
	phy$eventBranchSegs <- eventBranchSegs; 	
	phy$tipMu <- tipMu;
	phy$meanTipMu = meanTipMu;
	if (type == 'diversification') {
		phy$type <- 'diversification';
	} 
	else {
		phy$type <- 'trait';	
	}
 	
	# adds new class: 'bammdata'
	class(phy) <- 'bammdata';
	return(phy);
}
