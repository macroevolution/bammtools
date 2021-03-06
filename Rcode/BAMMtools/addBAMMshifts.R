addBAMMshifts = function(ephy, method, index, cex=1, pch=21, col=1, bg=2, multi=FALSE) {
	if (!'bammdata' %in% class(ephy)) stop("Object ephy must be of class bammdata");
	lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv);
	
	if (!multi){
		op = par(no.readonly = TRUE);
		par(lastPP$pp);		
	}

	
	shiftnodes = getShiftNodesFromIndex(ephy, index);
	isShift = ephy$eventData[[index]]$node %in% shiftnodes;
	times = ephy$eventData[[index]]$time[isShift];	
	
	if (method == 'phylogram') {
		if (max(lastPP$xx) <= 1) {
			XX = times/max(branching.times(as.phylo.bammdata(ephy)));
		} else {
			XX = times;
		}
		YY = lastPP$yy[shiftnodes];
	} else if (method == 'polar') {
		rb = lastPP$rb;
		XX = (rb+times/max(branching.times(as.phylo.bammdata(ephy)))) * cos(lastPP$theta[shiftnodes]);
		YY = (rb+times/max(branching.times(as.phylo.bammdata(ephy)))) * sin(lastPP$theta[shiftnodes]);		
	}	
	points(XX,YY,pch=pch,cex=cex,col=col,bg=bg);
	if (!multi) {
		par(op);		
	}

}