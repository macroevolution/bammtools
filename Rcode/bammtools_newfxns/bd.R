#
# Expected diversity of a time inhomogeneous birth death process
#
exd <- function(Tend, Tstart, args) {
	#args[1] t0
	#args[2] lam.init
	#args[3] lam.shape
	#args[4] mu.init
	#args[5] mu.shape
	rt <- function(Tend, Tstart, args) {
		Tmax <- Tend   - args[1];
		Tmin <- Tstart - args[1];
		if (abs(args[3]) > 0.00001)
			lam <- (args[2]/args[3])*(exp(args[3]*Tmax) - exp(args[3]*Tmin))
		else
			lam <- args[2]*(Tmax-Tmin);
		if (abs(args[5]) > 0.00001)
			mu  <- (args[4]/args[5])*(exp(args[5]*Tmax) - exp(args[5]*Tmin)) 
		else
			mu  <- args[4]*(Tmax-Tmin);
		return (mu - lam);	
	}
	pt <- function(Tend, Tstart, args) {
		integrand <- function(tau, Tstart, args) {
			(args[4]*exp(args[5]*(tau-args[1]))) * exp(rt(tau, Tstart, args));
		}
		res <- integrate(integrand, lower = Tstart, upper = Tend, Tstart=Tstart, args=args);
		if (res$message == "OK")
			return ( 1/(1 + res$value) )	
		else 
			stop("Integration failure");	
	} 
	b <- function(Tend, Tstart, args) {
		1 - pt(Tend, Tstart, args)*exp(rt(Tend, Tstart, args));
	}
	p <- 1-b(Tend, Tstart, args); # parameter of geometric distribution
	if (p == 0)
		return (0);
	return (1/p);                 # expected value of geometric distribution
}

#
# Calculate rate through time plot conditional on expected diversity
#
rtt <- function(ephy, nslices = 100, shift.include = NULL, shift.exclude = NULL) {
	if (is.null(shift.include))
		shift.include = c(ephy$edge[1,1], ephy$edge[,2]);
	if (!is.null(shift.exclude))
		shift.include = setdiff(shift.include, shift.exclude);
	if (ephy$type == "trait")
		stop("rtt is specific to diversification BAMM");
	root <- ephy$Nnode+2;
	bt <- branching.times(as.phylo.bammdata(ephy));
	tH <- max(bt);
	st <- min(tH-bt[as.character(intersect(shift.include,ephy$edge[,1]))]);
	tvec <- seq(st, tH, length.out = nslices);
	rt <- matrix(0, length(ephy$eventData), length(tvec));
	for (i in 1:length(ephy$eventData)) {
		ed <- ephy$eventData[[i]];
		ed <- ed[ed$node %in% shift.include,];
		if (nrow(ed) < 1)
			next;
		wts <- matrix(0,nrow(ed),ncol(rt));
		rts <- matrix(0,nrow(ed),ncol(rt));
		for (j in 1:nrow(ed)) {
			k <- which(tvec >= ed[j,2])[1];
			if (ed[j,1] == root)
				wts[j,k:ncol(rt)] <- 2*sapply(tvec[k:ncol(rt)], exd, ed[j,2], as.numeric(ed[j,2:6]))
			else 
				wts[j,k:ncol(rt)] <- sapply(tvec[k:ncol(rt)], exd, ed[j,2], as.numeric(ed[j,2:6]));			
			rts[j,k:ncol(rt)] <- sapply(tvec[k:ncol(rt)]-ed[j,2], exponentialRate, ed[j,3],ed[j,4]);
		}
		rt[i,] <- colSums(rts*wts)/colSums(wts);
	}
	return (list(rates=rt,times=tvec));
}

#
# how does it compare?
#

data(whales,events.whales);
ed <- getEventData(whales,events.whales,0.2);

r1 <- rtt(ed);
r2 <- getRateThroughTimeMatrix(ed);
	
plot(r1$times, colMeans(r1$rates), type='l', col=2, ylim = c(0.05,0.7), xlab="time", ylab="speciation rate");
lines(r2$times, colMeans(r2$beta), col=3);
legend('topleft',legend=c("analytic","posterior"),col=c(2,3),lty=c(1,1),bty="n")

