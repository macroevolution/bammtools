cohorts <- function(x, ephy, col, pal, lwd = 1, ofs = 0, use.plot.bammdata = TRUE, useraster = FALSE, LARGE = 500) {
	op <- par(no.readonly = TRUE);
	figs <- matrix(c(0,0.2,0.8,1,
	                 0.2,0.95,0.8+ofs,1,
	                 0,0.2-ofs,0,0.8,
	                 0.2,0.95,0,0.8,
	                 0.95,1,0.25,0.75
	                 ), byrow=TRUE,
	               nrow=5, ncol=4);
	if (dim(x)[1] > LARGE)
		useraster <- TRUE;
	if (missing(pal))
		pal <- "RdYlBu";
	if (missing(col))
		col <- colorRampPalette(get("palettes",.colorEnv)[["RdYlBu"]])(64);
	ncolors <- length(col);
	if (use.plot.bammdata) {               
		par(fig = figs[2,], new=FALSE, mar = c(0,0,1,4));
		plot(ephy, pal=pal,lwd=lwd,direction="downwards");
		par(fig = figs[3,], new=TRUE, mar = c(5,1,0,0));
		plot(ephy,pal=pal,lwd=lwd,direction="rightwards")
		par(fig = figs[4,], new=TRUE, mar = c(5,0,0,4));
		plot(0,0,type="n",axes=FALSE,ann=FALSE,xlim=c(0,1),ylim=c(0,1))
		image(x,axes=FALSE,xlab="",ylab="",col=col,xlim=c(0,1),ylim=c(0,1),add=TRUE,useRaster=useraster);
	}
	else {
		phy <- as.phylo.bammdata(ephy);
		bt <- max(branching.times(phy));
		par(fig = figs[2,], new=FALSE, mar = c(0,0,1,4));
		plot.phylo(phy,edge.width=lwd,direction="downwards",show.tip.label=FALSE,x.lim=c(1,length(phy$tip.label)),y.lim=c(0,bt));
		par(fig = figs[3,], new=TRUE, mar = c(5,1,0,0));
		plot.phylo(phy,edge.width=lwd,direction="rightwards",show.tip.label=FALSE,y.lim=c(1,length(phy$tip.label)),x.lim=c(0,bt));
		par(fig = figs[4,], new=TRUE, mar = c(5,0,0,4));
		gl <- 1:(length(ephy$tip.label)+1);
		plot(0,0,type="n",axes=FALSE,ann=FALSE,xlim=c(1,length(gl)-1),ylim=c(1,length(gl)-1))
		image(gl,gl,x,axes=FALSE,xlab="",ylab="",col=col,xlim=c(1,length(gl)-1),ylim=c(1,length(gl)-1),add=TRUE,useRaster=useraster);
	}
	barLegend(col, quantile(seq(min(x),max(x),length.out=ncolors+1),probs=seq(min(x),max(x),length.out=ncolors+1)),fig=figs[5,],side=2);
	par(op);
}
