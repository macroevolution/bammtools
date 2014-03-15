library(BAMMtools)
library(nlme)

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


res <- gls(y ~ x, data=dff, correlation=corBAMM(ephy=ed2))
# fails.






