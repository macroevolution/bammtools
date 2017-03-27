############################################
#	Internal function called by plot.bammdata(...)
#
#



colorMap <- function(x, pal, breaks, logcolor = FALSE, color.interval = NULL) {
    dpal <- get("palettes", envir = .colorEnv);
	NCOLORS <- length(breaks) - 1;
	if (length(pal) >= 3) {
		colpalette <- colorRampPalette(pal,space='Lab')(NCOLORS);	
	} else if (pal %in% names(dpal)) {
		colpalette <- colorRampPalette(dpal[[pal]],space='Lab')(NCOLORS);
	} else if (tolower(pal) == "temperature") {
		colpalette <- gplots::rich.colors(NCOLORS);	
	} else if (tolower(pal) == "terrain") {
		colpalette <- terrain.colors(NCOLORS);
	} else {
		stop("Unrecognized color palette specification");
	}
	
	if (!is.null(color.interval)) {
		if (length(color.interval) != 2) {
			stop("Color interval must have 2 values.");
		}
		if (is.na(color.interval[1])) {
			color.interval[1] <- min(breaks);
		}
		if (is.na(color.interval[2])) {
			color.interval[2] <- max(breaks);
		}
		
		#if color interval is contained within the range of supplied rates
		if (color.interval[1] >= min(breaks) & color.interval[2] <= max(breaks)) {
			goodbreaks <- intersect(which(breaks > color.interval[1]), which(breaks < color.interval[2]));
			topcolor <- colpalette[length(colpalette)];
			bottomcolor <- colpalette[1];
			
			NCOLORS <- length(goodbreaks) - 1;
			if (length(pal) >= 3) {
				colpalette2 <- colorRampPalette(pal,space='Lab')(NCOLORS);	
			} else if (pal %in% names(dpal)) {
				colpalette2 <- colorRampPalette(dpal[[pal]],space='Lab')(NCOLORS);
			} else if (tolower(pal) == "temperature") {
				colpalette2 <- gplots::rich.colors(NCOLORS);	
			} else if (tolower(pal) == "terrain") {
				colpalette2 <- terrain.colors(NCOLORS);
			}
			
			#replace colors in original color ramp
			colpalette[goodbreaks[1:(length(goodbreaks) - 1)]] <- colpalette2;
			colpalette[1:(goodbreaks[1] - 1)] <- bottomcolor;
			colpalette[goodbreaks[length(goodbreaks)]:length(colpalette)] <- topcolor;
		}

		#if color interval exceeds the range of color breaks
		if (color.interval[1] < min(breaks) | color.interval[2] > max(breaks)) {
			
			#generate new set of breaks for full range of rate values
			newbreaks <- seq(from = min(color.interval[1], min(breaks)), to = max(color.interval[2], max(breaks)), by = (breaks[2] - breaks[1]));
			newbreaks <- seq(from = min(color.interval[1], min(breaks)), to = max(color.interval[2], max(breaks)), length.out = length(newbreaks));
			
			# which breaks fall within the defined color.interval
			goodbreaks <- intersect(which(newbreaks >= color.interval[1]), which(newbreaks <= color.interval[2]));
			
			#generate colors for new breaks
			NCOLORS <- length(goodbreaks) - 1;
			if (length(pal) >= 3) {
				colpalette2 <- colorRampPalette(pal,space='Lab')(NCOLORS);	
			} else if (pal %in% names(dpal)) {
				colpalette2 <- colorRampPalette(dpal[[pal]],space='Lab')(NCOLORS);
			} else if (tolower(pal) == "temperature") {
				colpalette2 <- gplots::rich.colors(NCOLORS);	
			} else if (tolower(pal) == "terrain") {
				colpalette2 <- terrain.colors(NCOLORS);
			}
						
			# create new color palette that contains the color ramp 
			# within the color.interval
			# Fill in other slots with repeats of min or max color
			colpalette <- character(length(newbreaks) - 1);
			colpalette[goodbreaks[1:length(goodbreaks) - 1]] <- colpalette2;
			
			breaks <- newbreaks;
			
			if (any(colpalette == '')) {
				NAcol <- which(colpalette == '');
				nonNAcol <- which(colpalette != '');
				colFill <- sapply(NAcol, function(y) which.min(abs(y - nonNAcol)));
				colpalette[NAcol] <- colpalette[nonNAcol[colFill]];
			}
		}
	}
	
	kde <- density(x, from=min(x), to=max(x));
	colset <- numeric(length(x));
	coldens <- numeric(length(kde$x));
	for (i in 2:length(breaks)) {
        if (i == 2) {
            colset[x < breaks[2]] <- colpalette[1];
            coldens[kde$x < breaks[2]] <- colpalette[1];
        }
        else if (i == length(breaks)) {
            colset[x >= breaks[length(breaks)-1]] <- colpalette[length(breaks)-1];
            coldens[kde$x >= breaks[length(breaks)-1]] <- colpalette[length(breaks)-1];
        }
        else {
            colset[x >= breaks[i-1] & x < breaks[i]] <- colpalette[i-1];
        	coldens[kde$x >= breaks[i-1] & kde$x < breaks[i]] <- colpalette[i-1];
        }
    }
	coldens <- data.frame(kde$x,kde$y,coldens,stringsAsFactors=FALSE);	
	return(list(cols = colset, colsdensity = coldens, breaks = breaks, colpalette = colpalette));
}

# colorMap = function(x, pal, NCOLORS)
# {
	# dpal = c('BrBG','PiYG','PuOr','RdBu','RdGy','RdYlBu','RdYlGn','Spectral');
	# colset = numeric(length(x));
	# if(length(pal) == 3)
	# {
		# colpalette = colorRampPalette(pal,space='Lab')(NCOLORS);	
	# }
	# else if(pal %in% dpal)
	# {
		# colpalette = colorRampPalette(rev(brewer.pal(3,pal)),space='Lab')(NCOLORS);
	# }
	# else if(pal == 'temperature')
	# {
		# colpalette = gplots::rich.colors(NCOLORS);	
	# }
	
	# bks = quantile(x, seq(0,1,length.out=(NCOLORS+1)));
	# for(i in 2:length(bks))
	# {
		# if(i == 2)
		# {
			# colset[x < bks[2]] = colpalette[1];	
		# }
		# else if(i == length(bks))
		# {
			# colset[x >= bks[length(bks)-1]] = colpalette[length(bks)-1];
		# }
		# else
		# {
			# colset[x >= bks[i-1] & x < bks[i]] = colpalette[i-1];
		# }
	# }
	# return(colset);
# }
