library(BAMMtools)
library(nlme)
source('stable_corBAMM.R')

v <- read.tree('hackett_mcc.tre');
ed <- getEventData(v, 'eventdata_p50.txt', nsamples=200);

# Species-specific rates of evolution of reproductive isolation
ri <- read.csv('avian_speciesRIvals.txt', header=T);

goodsp <- intersect(v$tip.label, ri$sp)

# vector of species RI values
rivals <- ri$postri
names(rivals) <- ri$sp

# Get subtree from full avian tree with JUST taxa 
# in RI dataset
ed2 <- subtreeBAMM(ed, tips = names(rivals))

# Sorting..
rivals <- rivals[ed2$tip.label]

#########

# Now, extract tip speciation rates for same set of taxa
tiprates <- getTipRates(ed2)$lambda.avg
tiprates <-tiprates[names(rivals)]


# Build dataframe for GLS analysis
dff <- data.frame(x = rivals, y = tiprates)
rownames(dff) <- names(tiprates)




# test:
cm <- getCohortMatrix(ed2)

res <- gls(y ~ x, data=dff, correlation=corBAMM(ephy=ed2))

# 
res <- gls(y ~ x, data=dff, correlation=corPagel(0.5, as.phylo(ed2)))


### effective degrees of freedom:
cm <- vcv.phylo(as.phylo.bammdata(ed2), corr=T)
cm <- getCohortMatrix(ed2)

# Approach from Fraedrich et al. 1999
ee <- eigen(cm)
244^2 / (sum(ee$values^2))

 

# fudge to address the singularity of the bamm-correlation matrix:
cm[cm > 0.9999] <- 0.9999
diag(cm) <- rep(1, nrow(cm))
icm <- solve(cm) # matrix is computationally singular without the preceding fudge

# Estimator from Vallejos and Moreno:

ov <- t(rep(1, nrow(cm)))
ov %*% icm %*% t(ov)

###### Should be able to do this with simulations:

library(mvtnorm)


svec <- 1:244
REPS <- 10000

mm <- matrix(0, nrow=REPS, ncol=length(svec))
fx <- function(x) mean(rnorm(x))

for (i in 1:ncol(mm)){
	cat(i, '\n')
	mm[,i] <- sapply(rep(svec[i],REPS), fx)
	
}

vv <- apply(mm, 2, var)

## Now for observed:

ss <- rmvnorm(10000, sigma=cm)
v_obs <- rowMeans(ss)
var(v_obs);
max(which(var(v_obs) < vv))

## effective size of roughly ~12, similar to eigenvalue effective size.





