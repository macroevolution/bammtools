redirect <- function(coord, theta) {
	rot <- function(x, theta) {
		R <- matrix(c(cos(theta),-sin(theta),sin(theta),cos(theta)),byrow=TRUE,2,2);
		R%*%x;
	}
	tmp <- coord;
	tmp[,1:2] <- t(apply(coord[,1:2],1,rot,theta));
	tmp[,3:4] <- t(apply(coord[,3:4],1,rot,theta));
	return (tmp);
}

##' @title Plot \code{BAMM}-estimated macroevolutionary rates on a phylogeny
##'
##' @description \code{plot.bammdata} plots a phylogenetic tree from a
##'     \code{bammdata} object and colors each branch by the estimated rate of
##'     speciation, extinction, or trait evolution. Rates are not assumed to
##'     be constant in time, and the function can plot continuously-varying
##'     rates along individual branches.
##'
##' @param x An object of class \code{bammdata}.
##' @param tau A numeric indicating the grain size for the calculations. See
##'     documentation for \code{\link{dtRates}}.  	
##' @param method A character string indicating the method for plotting the
##'     phylogenetic tree. \code{method = "phylogram"} (default) plots the
##'     phylogenetic tree using rectangular coordinates.
##'     \code{method = "polar"} plots the phylogenetic tree using polar
##'     coordinates.
##' @param xlim A numeric vector of coordinates for the x-axis endpoints.
##'     Defaults to \code{NULL}, in which case they are calculated
##'     automatically. The x-axis corresponds to time when the phylogeny is
##'     plotted facing to the left or to the right. The time at the root
##'     equals zero.  
##' @param ylim A numeric vector of coordinates for the y-axis endpoints.
##'     Defaults to \code{NULL}, in which case they are calculated
##'     automatically. Tips are plotted at integer values beginning from zero
##'     and stepping by one when the phylogeny is plotted facing to the left
##'     or to the right.
##' @param vtheta A numeric indicating the angular separation (in degrees) of
##'     the first and last terminal nodes. Ignored if
##'     \code{method = "phylogram"}.
##' @param rbf A numeric indicating the length of the root branch as a
##'     fraction of total tree height. Ignored if \code{method = "phylogram"}.
##' @param show A logical indicating whether or not to plot the tree. Defaults
##'     to \code{TRUE}.
##' @param labels A logical indicating whether or not to plot the tip labels.
##'     Defaults to \code{FALSE}.
##' @param legend A logical indicating whether or not to plot a legend for
##'     interpreting the mapping of evolutionary rates to colors. Defaults to
##'     \code{FALSE}.
##' @param spex A character string indicating what type of macroevolutionary
##'     rates should be plotted. "s" (default) indicates speciation rates, "e"
##'     indicates extinction rates, and "netdiv" indicates net diversification
##'     rates. Ignored if \code{ephy$type = "trait"}.  	
##' @param lwd A numeric specifying the line width for branches.
##' @param cex A numeric specifying the size of tip labels.
##' @param pal A character string or vector of mode character that describes
##'     the color palette. See Details for explanation of options.
##' @param mask An optional integer vector of node numbers specifying branches
##'     that will be masked with \code{mask.color} when plotted.	
##' @param mask.color The color for the mask.  	
##' @param colorbreaks A numeric vector of percentiles delimiting the bins for
##'     mapping rates to colors. If \code{NULL} (default) bins are calculated
##'     from the rates that are passed with the \code{bammdata} object.  	
##' @param logcolor Logical. Should colors be plotted on a log scale.
##' @param breaksmethod Method used for determining color breaks. See help
##'     file for \code{\link{assignColorBreaks}}.
##' @param color.interval Min and max value for the mapping of rates. If
##'     \code{NULL}, then min and max are inferred from the data. NA can also
##'     be supplied for one of the two values. See details. 
##' @param JenksSubset If \code{breaksmethod = "jenks"}, the number of
##'     regularly spaced samples to subset from the full rates vector. Only
##'     relevant for large datasets. See help file for
##'     \code{\link{assignColorBreaks}}.
##' @param par.reset A logical indicating whether or not to reset the
##'     graphical parameters when the function exits. Defaults to
##'     \code{FALSE}. 
##' @param direction A character string. Options are "rightwards",
##'     "leftwards", "upwards", and "downwards", which determine the
##'     orientation of the tips when the phylogeny plotted.  	
##' @param \dots Further arguments passed to \code{par}.
##'
##' @details To calculate rates, each branch of the phylogeny is discretized
##'     into a number of small segments, and the mean of the marginal
##'     posterior density of the rate of speciation/extinction or trait
##'     evolution is calculated for each such segment. Rates are mapped to
##'     colors such that cool colors represent slow rates and warm colors
##'     represent fast rates. When the tree is plotted each of these small
##'     segments is plotted so that changes in rates through time and shifts
##'     in rates are visible as gradients of color. The \code{spex} argument
##'     determines the type of rate that will be calculated. \code{spex = "s"}
##'     will plot speciation rates, \code{spex = "e"} will plot extinction
##'     rates, and \code{spex = "netdiv"} will plot diversification rates
##'     (speciation - extinction). Note that if \code{x$type = "trait"} the
##'     \code{spex} argument is ignored and rates of phenotypic evolution are
##'     plotted instead. If \code{legend = TRUE} the function will plot a
##'     legend that contains the mapping of colors to numerical values.  
##'
##'     A number of color palettes come built in with \code{BAMMtools}.
##'     Color-blind friendly options include:
##'     \itemize{
##'         \item{BrBG}
##'         \item{PiYG}
##'         \item{PRGn}
##'         \item{PuOr}
##'         \item{RdBu}
##'         \item{RdYlBu}
##'         \item{BuOr}
##'         \item{BuOrRd}
##'         \item{DkRdBu}
##'         \item{BuDkOr}
##'         \item{GnPu}
##'     }
##'     Some color-blind unfriendly options include:
##'     \itemize{
##'         \item{RdYlGn}
##'         \item{Spectral}
##'         \item{temperature}
##'         \item{terrain}
##'     }
##'     Some grayscale options include:
##'     \itemize{
##'         \item{grayscale}
##'         \item{revgray}
##'     }
##'     For more information about these color palettes visit
##'     \url{https://colorbrewer2.org/} and
##'     \url{https://pjbartlein.github.io/datagraphics/color_scales.html} or
##'     use the help files of the R packages \code{RColorBrewer} and
##'     \code{dichromat}.
##'
##'     Additionally, any vector of valid named colors may also be used. The
##'     only restriction is that the length of this vector be greater than or
##'     equal to three (you can provide a single color, but in this case the
##'     entire tree will be assigned the same color). The colors should be
##'     ordered from cool to warm as the colors will be mapped from low rates
##'     to high rates in the order supplied (e.g. \code{pal=c("darkgreen",
##'     "yellow2", "red")}). The option \code{pal = "temperature"} uses the
##'     \code{rich.colors} function written by Arni Magnusson for the R
##'     package \code{gplots}.
##'
##'     Internally \code{plot.bammdata} checks whether or not rates have been
##'     calculated by looking for a component named "dtrates" in the
##'     \code{bammdata} object. If rates have not been calculated
##'     \code{plot.bammdata} calls \code{dtRates} with \code{tau}. Specifying
##'     smaller values for \code{tau} will result in smoother-looking rate
##'     changes on the tree. Note that smaller values of \code{tau} require
##'     more computation. If the \code{colorbreaks} argument
##'     is \code{NULL} a map of rates to colors is also made by calling
##'     \code{assignColorBreaks} with \code{NCOLORS = 64}. A user supplied
##'     \code{colorbreaks} argument can be passed as well. This allows one to
##'     plot parts of a tree while preserving the map of rates to colors that
##'     was made using rates for the entire tree.
##'
##'     If color.interval is defined, then those min and max values override
##'     the automatic detection of min and max. This might be useful if some
##'     small number of lineages have very high or very low rates, such that
##'     the map of colors is being skewed towards these extremes, resulting in
##'     other rate variation being drowned out. If specified, the color ramp
##'     will be built between these two color.interval values, and the rates
##'     outside of the color interval range will be set to the highest and
##'     lowest color. The total number of colors will also be increased such
##'     that 64 color bins are found within the color.interval.
##'
##'     If \code{plot.bammdata} is called repeatedly with the same
##'     \code{bammdata} object, computation can be reduced by first calling
##'     \code{dtRates} in the global environment.
##'
##' @return Returns (invisibly) a list with three components. 
##'     \itemize{
##'         \item{coords} {A matrix of plot coordinates. Rows correspond to
##'             branches. Columns 1-2 are starting (x,y) coordinates of each
##'             branch and columns 3-4 are ending (x,y) coordinates of each
##'             branch. If \code{method = "polar"} a fifth column gives the
##'             angle(in radians) of each branch.}
##'         \item{colorbreaks} {A vector of percentiles used to group
##'             macroevolutionary rates into color bins.}
##'         \item{colordens} {A matrix of the kernel density estimates (column
##'             2) of evolutionary rates (column 1) and the color (column 3)
##'             corresponding to each rate value.}
##'     }
##'
##' @source \url{https://colorbrewer2.org/},
##'     \url{https://pjbartlein.github.io/datagraphics/color_scales.html}
##'
##' @author Mike Grundler, Pascal Title
##'
##' @seealso \code{\link{dtRates}}, \code{\link{addBAMMshifts}},
##'     \code{\link{assignColorBreaks}}, \code{\link{subtreeBAMM}},
##'     \code{\link{colorRampPalette}}
##'
##' @examples
##' data(whales, events.whales)
##' ed <- getEventData(whales, events.whales, burnin=0.25, nsamples=500)
##'
##' # The first call to plot.bammdata
##' # No calculations or assignments of rates have been made
##' plot(ed, lwd = 3, spex = "s") # calls dtRates & assignColorBreaks
##'
##' # Compare the different color breaks methods
##' par(mfrow=c(1,3))
##' plot(ed, lwd = 3, spex = "s", breaksmethod = "linear")
##' title(main="linear")
##' plot(ed, lwd = 3, spex = "s", breaksmethod = "quantile")
##' title(main="quantile")
##' plot(ed, lwd = 3, spex = "s", breaksmethod = "jenks")
##' title(main="jenks")
##'
##' \dontrun{
##' # now plot.bammdata no longer calls dtRates
##' ed <- dtRates(ed, tau = 0.01)
##' xx <- plot(ed, lwd = 3, spex = "s")
##'
##' # you can plot subtrees while preserving the original 
##' # rates to colors map by passing the colorbreaks object as an argument
##' sed <- subtreeBAMM(ed, node = 103)
##' plot(sed, lwd = 3, colorbreaks = xx$colorbreaks)
##' sed <- subtreeBAMM(ed, node = 140)
##' plot(sed, lwd = 3, colorbreaks = xx$colorbreaks)
##' # note how if we do not pass colorbreaks the map is 
##' # no longer relative to the rest of the tree and the plot is quite
##' # distinct from the original
##' plot(sed, lwd = 3)
##'
##' # if you want to change the value of tau and the rates to colors map for
##' # the entire tree
##' ed <- dtRates(ed, tau = 0.002)
##' xx <- plot(ed, lwd = 3, spex = "s")
##' # now you can re-plot the subtrees using this finer tau partition
##' sed <- subtreeBAMM(ed, node = 103)
##' sed <- dtRates(sed, 0.002)
##' plot(sed, lwd = 3, colorbreaks = xx$colorbreaks)
##' sed <- subtreeBAMM(ed, node = 140)
##' sed <- dtRates(sed, 0.002)
##' plot(sed, lwd = 3, colorbreaks = xx$colorbreaks)
##'
##' # multi-panel plotting and adding shifts of specific posterior samples
##' par(mfrow=c(2,3))
##' samples <- sample(1:length(ed$eventData), 6)
##' ed <- dtRates(ed, 0.005)
##' # individual plots will have a color map relative to the mean
##' xx <- plot(ed, show=FALSE)
##' for (i in 1:6) {
##'     ed <- dtRates(ed, 0.005, samples[i])
##'     plot(ed, colorbreaks=xx$colorbreaks)
##'     addBAMMshifts(ed,index=samples[i],method="phylogram", par.reset=FALSE)	
##' }
##' dev.off()
##'
##' # color options
##' ed <- dtRates(ed,0.01)
##' plot(ed, pal="temperature",lwd=3)
##' plot(ed, pal="terrain",lwd=3)
##' plot(ed, pal=c("darkgreen","yellow2","red"),lwd=3)
##' plot(ed,method="polar",pal="Spectral", lwd=3)
##' plot(ed,method="polar",pal="RdYlBu", lwd=3)}

##' @keywords models graphics
##' @rdname plot
##' @aliases plot.bammdata
##' @export
##' @export plot.bammdata

plot.bammdata <- function (x, tau = 0.01, method = "phylogram", xlim = NULL, ylim = NULL, vtheta = 5, rbf = 0.001, show = TRUE, labels = FALSE, legend = FALSE, spex = "s", lwd = 1, cex = 1, pal = "RdYlBu", mask = integer(0), mask.color = gray(0.5), colorbreaks = NULL, logcolor = FALSE, breaksmethod = "linear", color.interval = NULL, JenksSubset = 20000, par.reset = FALSE, direction = "rightwards", ...) {
    if (inherits(x, "bammdata")) {
    	if (attributes(x)$order != "cladewise") {
    		stop("Function requires tree in 'cladewise' order");
    	}
        phy <- as.phylo.bammdata(x);
    }
    else stop("Object ephy must be of class bammdata");
    
    if (!spex %in% c('s','e','netdiv')) {
    	stop("spex must be 's', 'e' or 'netdiv'.");
    }
    
    if (length(pal) == 1 && !pal %in% names(get("palettes", envir=.colorEnv)) && pal != "temperature" && pal != "terrain")
    	pal <- rep(pal, 3)
    else if (length(pal) == 2)
    	pal <- c(pal, pal[2]);
    
    if (breaksmethod == 'linear' & !is.null(color.interval)) {
        if (length(color.interval) != 2) {
            stop("color.interval must be a vector of 2 numeric values.");
        }
    }
    
    if (!is.binary.phylo(phy)) {
        stop("Function requires fully bifurcating tree");
    }
    if (any(phy$edge.length == 0)) {
        warning("Tree contains zero length branches. Rates for these will be NA and coerced to zero");
    }
    if (!("dtrates" %in% names(x))) {
        x <- dtRates(x, tau);
    }
    
    NCOLORS <- 64;
	
	if (!is.null(color.interval)) {
    	# change the number of breaks such that the range of color.interval 
    	# is equivalent in terms of number of colors to the full range
    	# this way we preserve good resolution
        # Here we will ensure that NCOLORS bins occur within the color.interval
    	
   		if (x$type == "trait") {
    		ratesRange <- range(x$dtrates$rates);
    	} else if (x$type == "diversification") {
    		if (tolower(spex) == "s") {
    			ratesRange <- range(x$dtrates$rates[[1]]);
    		} else if (tolower(spex) == "e") {
    			ratesRange <- range(x$dtrates$rates[[2]]);
    		} else if (tolower(spex) == "netdiv") {
    			ratesRange <- range(x$dtrates$rates[[1]] - x$dtrates$rates[[2]]);
    		}
    	}	
    	if (all(!is.na(color.interval))) {
    		brks <- seq(min(color.interval[1], ratesRange[1]), max(color.interval[2], ratesRange[2]), length.out = (NCOLORS+1));
    		intervalLength <- length(which.min(abs(color.interval[1] - brks)) : which.min(abs(color.interval[2] - brks)));
    	} else if (is.na(color.interval[1])) {
    		brks <- seq(ratesRange[1], max(color.interval[2], ratesRange[2]), length.out = (NCOLORS+1));
    		intervalLength <- length(1 : which.min(abs(color.interval[2] - brks)));
    	} else if (is.na(color.interval[2])) {
    		brks <- seq(min(color.interval[1], ratesRange[1]), ratesRange[2], length.out = (NCOLORS+1));
    		intervalLength <- length(which.min(abs(color.interval[1] - brks)) : length(brks));
    	}
    	NCOLORS <- round((NCOLORS ^ 2) / intervalLength)    	
    }

    if (is.null(colorbreaks)) {
   	    colorbreaks <- assignColorBreaks(x$dtrates$rates, NCOLORS, spex, logcolor, breaksmethod, JenksSubset);
    }
    if (x$type == "trait") {
    	colorobj <- colorMap(x$dtrates$rates, pal, colorbreaks, logcolor, color.interval);
    }
    else if (x$type == "diversification") {
        if (tolower(spex) == "s") {
            colorobj <- colorMap(x$dtrates$rates[[1]], pal, colorbreaks, logcolor, color.interval);
        }
        else if (tolower(spex) == "e") {
            colorobj <- colorMap(x$dtrates$rates[[2]], pal, colorbreaks, logcolor, color.interval);
        }
        else if (tolower(spex) == "netdiv") {
            colorobj <- colorMap(x$dtrates$rates[[1]] - x$dtrates$rates[[2]], pal, colorbreaks, logcolor, color.interval);
        }
    }
    else {
   	    stop("Unrecognized/corrupt bammdata class. Type does not equal 'trait' or 'diversification'");	
    }
    edge.color <- colorobj$cols;
#    if (is.ultrametric(phy))    
#    	tH <- max(branching.times(phy))
#    else
#    	tH <- max(NU.branching.times(phy));
	tH <- max(x$end);
    phy$begin <- x$begin;
    phy$end <- x$end;
    tau <- x$dtrates$tau;
    if (method == "polar") {
        ret <- setPolarTreeCoords(phy, vtheta, rbf);
        rb <- tH * rbf;
        p <- mkdtsegsPolar(ret$segs[-1,], tau, x$edge);
    }
    else if (method == "phylogram") {
        ret <- setPhyloTreeCoords(phy);
        p <- mkdtsegsPhylo(ret$segs[-1,], tau, x$edge);
    }
    else {
        stop("Unimplemented method");
    }
    x0 <- c(ret$segs[1,1], p[, 1]);
    x1 <- c(ret$segs[1,3], p[, 2]);
    y0 <- c(ret$segs[1,2], p[, 3]);
    y1 <- c(ret$segs[1,4], p[, 4]);
    offset <- table(p[, 5])[as.character(unique(p[, 5]))];
    if (length(mask)) {
   	    edge.color[p[,5] %in% mask] <- mask.color;
    }
    arc.color <- c(edge.color[1], edge.color[match(unique(p[, 5]), p[, 5]) + offset]);
    edge.color <- c(edge.color[1], edge.color);
    if (show) {
    	    op <- par(no.readonly = TRUE);
        if (length(list(...))) {
            par(...);
        }
        if (legend) {
            #par(fig=c(0,0.9,0,1));
            par(mar = c(5, 4, 4, 5))
        }
        plot.new();
        ofs <- 0;
        if (labels) {
        	if (method == "phylogram")
	            ofs <- max(nchar(phy$tip.label) * 0.03 * cex * tH)
        	else
        		ofs <- max(nchar(phy$tip.label) * 0.03 * cex);
        }
        if (method == "polar") {
            if (is.null(xlim) || is.null(ylim)) {
                if (is.null(xlim))
                    xlim = c(-1, 1) + c(-rb, rb) + c(-ofs, ofs)
                if (is.null(ylim))
                    ylim = c(-1, 1) + c(-rb, rb) + c(-ofs, ofs) 
            }
            plot.window(xlim = xlim, ylim = ylim, asp = 1);
            segments(x0, y0, x1, y1, col = edge.color, lwd = lwd, lend = 2);
            arc(0, 0, ret$arcs[, 1], ret$arcs[, 2], c(rb, rb + phy$end/tH), border = arc.color, lwd = lwd);
            if (labels) {
                for (k in 1:length(phy$tip.label)) {
                  text(ret$segs[-1, ][phy$edge[, 2] == k, 3],ret$segs[-1, ][phy$edge[, 2] == k, 4], phy$tip.label[k],cex = cex, srt = (180/pi) * ret$arcs[-1,][phy$edge[, 2] == k, 1], adj = c(0, NA));
                }
            }
        }
        if (method == "phylogram") {
            direction <- match.arg(direction, c("rightwards","leftwards","downwards","upwards"));
        	if (direction == "rightwards") {
        		bars <- redirect(cbind(x0,y0,x1,y1),0);
            	arcs <- redirect(ret$arcs,0);
            	bars[,c(1,3)] <- tH * bars[,c(1,3)];
            	arcs[,c(1,3)] <- tH * arcs[,c(1,3)];
            	
            	# xlim <- c(0, 1 + ofs);
        		# ylim <- c(0, phy$Nnode * 1/(phy$Nnode + 1));
            	
            	ret$segs[-1, c(1,3)] <- tH * ret$segs[-1, c(1,3)]; 
            	        	
        	}
        	else if (direction == "leftwards") {
        		bars <- redirect(cbind(x0,y0,x1,y1),pi);
        		bars[,c(2,4)] <- abs(bars[,c(2,4)]);
            	arcs <- redirect(ret$arcs,pi);
            	arcs[,c(2,4)] <- abs(arcs[,c(2,4)]);
            	

            	bars[,c(1,3)] <- tH * bars[,c(1,3)];
            	arcs[,c(1,3)] <- tH * arcs[,c(1,3)];
            	
            	
            	ret$segs[-1, c(1,3)] <- -tH * ret$segs[-1, c(1,3)];
            	
				# xlim <- rev(-1*c(0, 1 + ofs));
				# ylim <- c(0, phy$Nnode * 1/(phy$Nnode + 1));
        	}
        	else if (direction == "downwards") {
        		bars <- redirect(cbind(x0,y0,x1,y1),-pi/2);
            	arcs <- redirect(ret$arcs,-pi/2);
            	
            	bars[,c(2,4)] <- tH * bars[,c(2,4)];
            	arcs[,c(2,4)] <- tH * arcs[,c(2,4)];
            	
            	
            	ret$segs <- redirect(ret$segs, -pi/2);
            	ret$segs[,c(2,4)] <- tH * ret$segs[,c(2,4)];
            	
            	# xlim <- c(0, phy$Nnode * 1/(phy$Nnode + 1));
            	# ylim <- rev(-1*c(0, 1 + ofs));	
        	}
        	else if (direction == "upwards") {
        		bars <- redirect(cbind(x0,y0,x1,y1),pi/2);
        		bars[,c(1,3)] <- abs(bars[,c(1,3)]);
            	arcs <- redirect(ret$arcs,pi/2);
            	arcs[,c(1,3)] <- abs(arcs[,c(1,3)]);
            	
            	bars[,c(2,4)] <- tH * bars[,c(2,4)];
            	arcs[,c(2,4)] <- tH * arcs[,c(2,4)];
            	
            	ret$segs <- redirect(ret$segs, pi/2);
            	ret$segs[,c(1,3)] <- abs(ret$segs[,c(1,3)]);
            	ret$segs[,c(2,4)] <- tH * ret$segs[,c(2,4)];
            	
        		# xlim <- c(0, phy$Nnode * 1/(phy$Nnode + 1));
        		# ylim <- c(0, 1 + ofs);
        	}
        	if (is.null(xlim) && direction == "rightwards") xlim <- c(0, tH + ofs);
        	if (is.null(xlim) && direction == "leftwards") xlim <- c(-(tH + ofs), 0);
        	if (is.null(ylim) && (direction == "rightwards" || direction == "leftwards")) ylim <- c(0, phy$Nnode);  
        	
        	if (is.null(xlim) && (direction == "upwards" || direction == "downwards")) xlim <- c(0, phy$Nnode);
        	if (is.null(ylim) && direction == "upwards") ylim <- c(0, tH + ofs);
        	if (is.null(ylim) && direction == "downwards") ylim <- c(-(tH + ofs), 0);  
        	
        	   
            plot.window(xlim = xlim, ylim = ylim);
            segments(bars[-1,1], bars[-1,2], bars[-1,3], bars[-1,4], col = edge.color[-1], lwd = lwd, lend = 2);
            isTip <- phy$edge[, 2] <= phy$Nnode + 1;
            isTip <- c(FALSE, isTip);
            segments(arcs[!isTip, 1], arcs[!isTip, 2], arcs[!isTip, 3], arcs[!isTip, 4], col = arc.color[!isTip], lwd = lwd, lend = 2);  
            if (labels) {
                if (direction == "rightwards")
 	                text(ret$segs[isTip, 3], ret$segs[isTip, 4], phy$tip.label[phy$edge[isTip[-1],2]], cex = cex, pos = 4, offset = 0.25)
                else if (direction == "leftwards")
                    text(ret$segs[isTip, 3], ret$segs[isTip, 4], phy$tip.label[phy$edge[isTip[-1],2]], cex = cex, pos = 2, offset = 0.25)
                else if (direction == "upwards")
                    text(ret$segs[isTip, 3], ret$segs[isTip, 4], phy$tip.label[phy$edge[isTip[-1],2]], cex = cex, pos = 4, srt = 90, offset = 0)
                else if (direction == "downwards")
                    text(ret$segs[isTip, 3], ret$segs[isTip, 4], phy$tip.label[phy$edge[isTip[-1],2]], cex = cex, pos = 2, srt = 90, offset = 0);
            }
        }
        # if (legend) {
            # #rateLegend(colorobj$colsdensity, logcolor);
            # if (is.null(color.interval)) {
            	# barLegend(pal, colorbreaks, fig=c(0.9,1,0.25,0.75), side=2);
        	# } else {
        		# barLegend(pal, colorbreaks, fig=c(0.9,1,0.25,0.75), side=2, colpalette=colorobj$colpalette);
        	# }
        # }
    }
    index <- order(as.numeric(rownames(ret$segs)));
    if (show) {
    	if (method == "phylogram") {
        	assign("last_plot.phylo", list(type = "phylogram", direction = direction, Ntip = phy$Nnode + 1, Nnode = phy$Nnode, edge = phy$edge, xx = ret$segs[index, 3], yy = ret$segs[index, 4], pp = par(no.readonly = TRUE)), envir = .PlotPhyloEnv);
		} else if (method == "polar") {
        	assign("last_plot.phylo", list(type = "fan", Ntip = phy$Nnode + 1, Nnode = phy$Nnode, edge = phy$edge, xx = ret$segs[index, 3], yy = ret$segs[index, 4], theta = ret$segs[index, 5], rb = rb, pp = par(no.readonly = TRUE)), envir = .PlotPhyloEnv);
		}
	if (legend) {
		addBAMMlegend(x = list(coords = ret$segs[-1, ], colorbreaks = colorobj$breaks, palette = colorobj$colpalette, colordens = colorobj$colsdensity), location = 'right')
	}
	}
    if (par.reset) {
        par(op);
    }
    invisible(list(coords = ret$segs[-1, ], colorbreaks = colorobj$breaks, palette = colorobj$colpalette, colordens = colorobj$colsdensity));
}










