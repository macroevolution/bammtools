cohorts <- function(x, ephy, col, pal, use.plot.bammdata = TRUE) {
	op <- par(no.readonly = TRUE);
	figs <- matrix(c(0,0.2,0.8,1,
	                 0.2,1,0.8,1,
	                 0,0.2,0,0.8,
	                 0.2,1,0,0.8
	                 ), byrow=TRUE,
	               nrow=4, ncol=4);
	if (missing(pal))
		pal <- "RdYlBu";
	if (missing(col))
		col <- colorRampPalette(get("palettes",.colorEnv)[["RdYlBu"]])(64);
	ncolors <- length(col);
	if (use.plot.bammdata) {               
		par(fig = figs[2,], new=FALSE, mar = c(0,0,1,4));
		plot(ephy, pal=pal,direction="downwards");
		par(fig = figs[3,], new=TRUE, mar = c(5,1,0,0));
		plot(ephy,pal=pal,direction="rightwards")
		par(fig = figs[4,], new=TRUE, mar = c(5,0,0,4));
		plot(0,0,type="n",axes=FALSE,ann=FALSE,xlim=c(0,1),ylim=c(0,1))
		image(x,axes=FALSE,xlab="",ylab="",col=col,xlim=c(0,1),ylim=c(0,1),add=TRUE);
	}
	else {
		phy <- as.phylo.bammdata(ephy);
		par(fig = figs[2,], new=FALSE, mar = c(0,0,1,4));
		plot.phylo(phy,direction="downwards",show.tip.label=FALSE,x.lim=c(1,length(phy$tip.label)));
		par(fig = figs[3,], new=TRUE, mar = c(5,1,0,0));
		plot.phylo(phy,direction="rightwards",show.tip.label=FALSE,y.lim=c(1,length(phy$tip.label)));
		par(fig = figs[4,], new=TRUE, mar = c(5,0,0,4));
		gl <- 1:(length(ephy$tip.label)+1);
		plot(0,0,type="n",axes=FALSE,ann=FALSE,xlim=c(1,length(gl)-1),ylim=c(1,length(gl)-1))
		image(gl,gl,x,axes=FALSE,xlab="",ylab="",col=col,xlim=c(1,length(gl)-1),ylim=c(1,length(gl)-1),add=TRUE);
	}
	barLegend(col, quantile(seq(min(x),max(x),length.out=ncolors+1),probs=seq(min(x),max(x),length.out=ncolors+1)));
	par(op);
}

data(whales,events.whales);
ed <- getEventData(whales,events.whales,0.1);
cormat <- getCohortMatrix(ed);

cohorts(cormat,ed);
cohorts(cormat,ed,use.plot.bammdata=FALSE);
cohorts(cormat,ed,col=terrain.colors(64));
cohorts(cormat,ed,pal="RdYlGn");
