## show.all.nodes: should plot core and "floater" nodes with different colors 
#	that can be specified by the user.
#    argument to specify layout , e.g., c(2,2) would specify a 2x2 panel?
#	 maybe print the number of shift configurations above the threshold that were 
#	 not plotted along with their frequency?
#	
#	add.freq.text = argument to add text to each plot indicating the frequency of the 
#		shift configuration

##' @title Plot distinct rate shift configurations on a phylogeny
##'
##' @description Plots a random distinct rate shift configuration sampled by
##'     \code{BAMM} on a phylogeny.
##'
##' @param x An object of class \code{bammshifts}.
##' @param ephy An object of class \code{bammdata}.
##' @param method A character string for which plotting method to use.
##'     "phylogram" uses rectangular coordinates. "polar" uses polar
##'     coordinates.
##' @param pal The color palette to use in \code{plot.bammdata}.
##' @param rank The rank of the core shift configuration to plot. For the
##'     default (\code{NULL}) a random configuration is chosen. 
##' @param index The posterior sample to plot. For the default (\code{NULL})
##'     a random sample is chosen.
##' @param spex A character string indicating what type of macroevolutionary
##'     rates should be plotted. "s" (default) indicates speciation rates, "e"
##'     indicates extinction rates, and 'netdiv' indicates net diversification
##'     rates. Ignored if \code{ephy$type = "trait"}.
##' @param legend Logical indicating whether to plot a legend.
##' @param add.freq.text A logical indicating whether the frequency of each
##'     sampled shift configuration should be added to each plot.
##' @param logcolor Logical. Should colors be plotted on a log scale.
##' @param breaksmethod Method used for determining color breaks. See help
##'     file for \code{\link{assignColorBreaks}}.
##' @param color.interval Min and max value for the mapping of rates. One of
##'     the two values can be \code{NA}. See details in
##'     \code{\link{plot.bammdata}} for further details. 
##' @param JenksSubset If \code{breaksmethod = "jenks"}, the number of
##'     regularly-spaced samples to subset from the full rates vector. Only
##'     relevant for large datasets. See help file for
##'     \code{\link{assignColorBreaks}}.
##' @param \dots Other arguments to \code{plot.bammdata}.
##'
##' @details A rate shift configuration is the set of nodes of the phylogeny
##'     where a shift occurs in the macroevolutionary rate dynamic of
##'     diversification or trait evolution. Each posterior sample is a
##'     potentially distinct rate shift configuration. Different
##'     configurations may imply different macroevolutionary scenarios. This
##'     function helps visualize the different distinct rate shift
##'     configurations sampled by \code{BAMM}.
##'
##'     A core shift configuration is defined by a set of nodes that have
##'     shift probabilities that are substantially elevated relative to what
##'     you expect under the prior alone. These core configurations are
##'     inferred in \code{\link{distinctShiftConfigurations}}. It is almost
##'     certain that more than one core shift configuration will be sampled by
##'     \code{BAMM}. Moreover, each core shift configuration may contain many
##'     subconfigurations. A subconfiguration contains the core shift node
##'     configuration and zero or more additional shift nodes that occur with
##'     low marginal probability.  
##'
##'     Points are added to the branches subtending the nodes of each rate
##'     configuration. The size of the point is proportional to the marginal
##'     probability that a shift occurs on a specific branch. If the
##'     instantaneous rate at a shift's origin represents an initial increase
##'     above the ancestral instantaneous rate the point is red. If the
##'     instantaneous rate at a shift's origin represents an initial decrease
##'     below the ancestral instantaneous rate the point is blue.
##'
##' @author Mike Grundler, Dan Rabosky
##'
##' @seealso \code{\link{distinctShiftConfigurations}},
##'     \code{\link{plot.bammdata}}
##'
##' @references \url{http://bamm-project.org}
##'
##' @examples
##' data(whales, events.whales)
##' 
##' ed <- getEventData(whales, events.whales, burnin=0.25, nsamples=500)
##' 
##' sc <- distinctShiftConfigurations(ed, expectedNumberOfShifts = 1,
##'                                   threshold = 5)
##' 
##' plot(sc, ed)
##' @export
plot.bammshifts <- function(x, ephy, method="phylogram", pal="RdYlBu", 
rank=NULL, index=NULL, spex="s", legend=TRUE, add.freq.text=TRUE, logcolor=FALSE, breaksmethod="linear", color.interval=NULL, JenksSubset=20000, ...) 
{
	if (class(x) != "bammshifts") {
		stop("arg sc must be of class 'bammshifts'");
	}        	
	if (class(ephy) != "bammdata") {
		stop("arg ephy must be of class 'bammdata'");
	}
	if (!spex %in% c('s', 'e', 'netdiv')) {
		stop("arg spex must be 's', 'e' or 'netdiv'.")
	}
	if ((spex == "e" || spex == "netdiv") && ephy$type == "trait") {
		warning("arg spex not meaningful for BAMMtrait");
		spex <- "s";
	}
	if (is.null(rank) && is.null(index)) {
		rank <- sample.int(length(x$shifts),1);
		index <- x$samplesets[[rank]][sample.int(length(x$samplesets[[rank]]),1)];
	}
	else if (!is.null(rank) && is.null(index)) {
		index <- x$samplesets[[rank]][sample.int(length(x$samplesets[[rank]]),1)];
	}
	else if (is.null(rank) && !is.null(index)) {
		rank <- {
			for (i in 1:length(x$sampleset)) {
			    if (index %in% x$sampleset[[i]]) {
				    break;
			    }
			}
			i;
		}
	}
	else {
		if (index > length(x$samplesets[[rank]])) {
			warning("arg index is not relative to the set of posterior samples of the given core shift configuration");
			index <- x$samplesets[[rank]][1];
		}
		else {
    		index <- x$samplesets[[rank]][index];
		}
	}
	
	par.reset <- TRUE;
	if (legend) {
		par.reset <- FALSE;
		m <- matrix(c(1,0,1,2,1,0),byrow=TRUE,nrow=3,ncol=2);
    	layout(m, widths=c(1,0.25));
		par(mar=c(7.1,1,7.1,1));
	}
	
	ephy <- dtRates(ephy,0.01);
	colorbreaks <- assignColorBreaks(ephy$dtrates$rates, spex=spex, logcolor=logcolor, method=breaksmethod, JenksSubset=JenksSubset);
	sed <- subsetEventData(ephy, index);
	plot.bammdata(sed, method=method, pal=pal, spex=spex, colorbreaks=colorbreaks, par.reset=par.reset, ...);

	shiftnodes <- getShiftNodesFromIndex(ephy, index);
	shiftnode_parents <- ephy$edge[match(shiftnodes,ephy$edge[,2],nomatch=0), 1];
    root <- (shiftnode_parents == (ephy$Nnode + 2));
    if (sum(root) > 0) {
        isShiftNodeParent <- integer(length(shiftnodes));
        isShiftNodeParent[root] <- 1;
        isShiftNodeParent[!root] <- sed$eventVectors[[1]][match(shiftnode_parents[!root], ephy$edge[,2])];
    } 
    else {
        isShiftNodeParent <- sed$eventVectors[[1]][match(shiftnode_parents, ephy$edge[,2])];
    }
    isShiftNode <- match(shiftnodes, sed$eventData[[1]]$node);
    time <- sed$eventData[[1]][isShiftNode, 2] - sed$eventData[[1]][isShiftNodeParent, 2];
    if (spex == "s") {
    	lam1 <- sed$eventData[[1]][isShiftNodeParent, 3]; 
    	lam2 <- sed$eventData[[1]][isShiftNodeParent, 4];
        AcDc <- exponentialRate(time, lam1, lam2) > sed$eventData[[1]][isShiftNode, 3];	
    }
    else if (spex == "e") {
    	mu1 <- sed$eventData[[1]][isShiftNodeParent, 5]; 
    	mu2 <- sed$eventData[[1]][isShiftNodeParent, 6];
    	AcDc <- exponentialRate(time, mu1, mu2) > sed$eventData[[1]][isShiftNode, 5];
    }
    else if (spex == 'netdiv') {
    	lam1 <- sed$eventData[[1]][isShiftNodeParent, 3]; 
    	lam2 <- sed$eventData[[1]][isShiftNodeParent, 4];
    	mu1 <- sed$eventData[[1]][isShiftNodeParent, 5]; 
    	mu2 <- sed$eventData[[1]][isShiftNodeParent, 6];
    	AcDc <- (exponentialRate(time, lam1, lam2)-exponentialRate(time, mu1, mu2)) > (sed$eventData[[1]][isShiftNode, 3]-sed$eventData[[1]][isShiftNode, 5]);
    }
    bg <- rep("blue", length(AcDc));
    bg[which(AcDc == FALSE)] <- "red";
	cex <- 0.75 + 5 * x$marg.probs[as.character(getShiftNodesFromIndex(ephy, index))];
	addBAMMshifts(sed, 1, method, cex=cex, bg=transparentColor(bg, 0.5),par.reset=par.reset);
	if (add.freq.text) {
		mtext(sprintf("core shift configuration: rank %i of %i", rank, length(x$shifts)),3,line=0);
		mtext(sprintf("sampled with frequency f = %.2g",x$frequency[rank]),3,line=-1.25);
		mtext(sprintf("showing subconfiguration %i of %i with this rank", match(index,x$samplesets[[rank]]), length(x$samplesets[[rank]])),1,line=-1.25);
	}
	if (legend) {
		par(mar=c(5.1,1.1,2.1,2.1));
		plot.new();			
		plot.window(xlim=c(0,1),ylim=c(0,1.25));
		points(rep(0.5,8),rev(seq(0,1.125,length.out=8)),pch=21,bg="white",cex = rev(0.75 + 5*seq(0.125,1,0.125)));
		mtext("marginal shift probability",cex=0.75,line=0);
		text(x=rep(0.8,8),y=rev(seq(0,1.125,length.out=8)),labels=sprintf("%10.3f",rev(seq(0.125,1,0.125))));	
	}
}


# plot.bammshifts <- function(sc, ephy, plotmax=9, method='phylogram', pal = 'RdYlBu', spex = "s", add.arrows = TRUE, add.freq.text = TRUE, use.plot.bammdata = TRUE, send2pdf = FALSE, ...){
	# if (class(sc) != 'bammshifts') {
		# stop('arg sc must be of class "bammshifts"');
	# }
	# if (class(ephy) != 'bammdata') {
		# stop('arg ephy must be of class "bammdata"');
	# }
	# if (plotmax > 9 && send2pdf == FALSE) {
	    # plotmax = 9;
	    # cat("arg plotmax coerced to 9\n");
	# }
	# mm <- min(c(length(sc$frequency), plotmax));
	# if (send2pdf) {
	    # pdf("shiftconfig.pdf");
	# }
	# if (mm == 1) {
	    # par(mfrow=c(1,1));
	# } else if (mm <= 2) {
		# par(mfrow=c(1,2));
	# } else if (mm <= 4) {
		# par(mfrow = c(2,2));
	# } else if (mm <= 6) {
		# par(mfrow=c(2,3));
	# } else {
		# par(mfrow=c(3,3));
	# }
	# cat("Omitted", max(length(sc$frequency),mm) - min(length(sc$frequency),mm), "plots\n");
	# if (use.plot.bammdata) {
    	# ephy = dtRates(ephy, 0.01);
	    # colorbreaks = assignColorBreaks(ephy$dtrates$rates,spex=spex);
	# }
	# tH = max(branching.times(as.phylo.bammdata(ephy)));
	# for (i in 1:mm) {
		# tmp <- subsetEventData(ephy, index=sc$sampleset[[i]]);
		# par(mar = c(2,2,2,2));
		# if (use.plot.bammdata) {
    		# plot.bammdata(tmp, method, pal=pal, colorbreaks=colorbreaks, ...);
		# }
		# else {
		    # if (method=="polar") method = "fan";
		    # plot.phylo(as.phylo.bammdata(ephy),type=method,show.tip.label=FALSE);
		# }
		# if (add.freq.text) {
			# mtext(paste("f =",signif(sc$frequency[i],2)),3);
			# mtext(paste(length(sc$sampleset[[i]]),"subconfigurations"),1,1.5);
		# }
		# box();
		# cex = 2 + 8 * sc$marg.probs[as.character(getShiftNodesFromIndex(ephy,sc$sampleset[[i]][1]))];
		# shiftnodes = getShiftNodesFromIndex(ephy,sc$sampleset[[i]][1]);
		# shiftnode_parents = ephy$edge[which(ephy$edge[,2] %in% shiftnodes),1];
		
		# isShiftNodeParent = integer(length(shiftnodes));		
		# root = (shiftnode_parents == ephy$Nnode+2);
		# if (sum(root) > 0) {
			# isShiftNodeParent[root] = 1;
		# }
		# if (sum(!root) > 0) {
			# isShiftNodeParent[!root] = tmp$eventVectors[[1]][which(ephy$edge[,2] %in% shiftnode_parents[!root])];
		# }	
		# isShiftNode = which(tmp$eventData[[1]]$node %in% shiftnodes);
		
		# AcDc = exponentialRate(tmp$eventData[[1]][isShiftNode, 2]-tmp$eventData[[1]][isShiftNodeParent, 2],tmp$eventData[[1]][isShiftNodeParent,3],tmp$eventData[[1]][isShiftNodeParent,4]) > tmp$eventData[[1]][isShiftNode, 3];
		
		# bg = rep("blue", length(AcDc));
		# bg[which(AcDc == FALSE)] = "red";
		
		# addBAMMshifts(tmp, method, 1, cex = cex, bg = transparentColor(bg,0.5), multi=TRUE);
		# if (add.arrows) {
			# lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv);
    	    # XX = lastPP$xx[sc$shifts[[i]]];
        	# YY = lastPP$yy[sc$shifts[[i]]];
        	# arrows(XX + 0.15*cos(3*pi/4),YY + 0.15*sin(3*pi/4), XX, YY, length=0.1,lwd=2);
        	# #points(XX, YY, pch = 21, cex = cex, col = 1, bg = transparentColor(bg,0.5)); 		
	    # }
	# }
	# if (send2pdf) dev.off();
# }

