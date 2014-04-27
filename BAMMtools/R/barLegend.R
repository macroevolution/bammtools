barLegend <- function(pal, colorbreaks) {
	pal <- colorRampPalette(get("palettes",envir=.colorEnv)[[pal]])(length(colorbreaks)-1);
	n <- length(pal);
	x <- seq(0,n,1)/n;
	x <- rep(x,each=2);
	x <- x[-c(1,length(x))];
	x <- matrix(x,ncol=2,byrow=TRUE);
	par(fig=c(0.9,1,0.25,0.75),new=TRUE,mar=c(0,0,0,0));
	plot.new();
	plot.window(ylim=c(0,1),xlim=c(-0.5,0.5));
	segments(y0=x[,1],x0=-0.5,y1=x[,2],x1=-0.5,col=pal,lwd=8,lend=2);
	tx <- colorbreaks[c("0%","50%","100%")];
	axis(4,at=c(0,0.5,1),labels=signif(tx,1),line=-3, las=1);
}
