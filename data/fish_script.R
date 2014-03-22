
library(BAMMtools);

edstring <- 'fish_events2.txt';
treestring <- 'fishtreeFinalPL_Feb8.tre';

v <- read.tree(treestring);
ed <- getEventData(v, edstring, burnin=0.1)


#colscheme <- plot.bammdata(ed, pal='temperature', lwd=0.5, legend=T)
acanth <- extract.clade(v, node = 7923);


# first, super slow - maybe this can speed up...
# subtreeBAMM has problems with dataset of this size
ed_acanth <- subtreeBAMM(ed, tips=acanth$tip.label);


#### I want to plot acanthomorph and non-acanthomorph fishes on common timescale.
# First, verify that rates are OK for full event data object:

rmat <- getRateThroughTimeMatrix(ed)
plotRateThroughTime(rmat)

### Non-acanthomorph fishes:
plotRateThroughTime(ed, node = 7923, nodetype = 'exclude', start.time = 225, plot=T);
# Plot does not start at 225 million years before present, but I think it should based on
# these arguments

# Acanthomorphs: fails due to missing values
plotRateThroughTime(ed, node = 7923, plot=T, start.time = 225);

# yes, this overshoots the crown age of the clade, but I think
# we want the function to work. 
#	It should should make the plot but just not plot where there are NA values.



# Trying something different. OK, this gives right timescale
#	... but line does weird thing at the end
# 	this will confuse users.
plotRateThroughTime(ed, node = 7923, plot=T, xlim=c(225, 0));


# Try something different yet again
rmat_A <- getRateThroughTimeMatrix(ed, node = 7923)
plotRateThroughTime(rmat_A, xlim=c(300,0))
# works, but again gives weird dip at end...

# Can we avoid dip?
rmat_A <- getRateThroughTimeMatrix(ed, node = 7923, end.time = 1)

# no, this now gives an NA-based error message.
plotRateThroughTime(rmat_A, xlim=c(300,0))




# Other issues
#
# start.time argument is inconsistent 
#		between plotRateThroughTime and getRateThroughTimeMatrix
#
# Should probably make this argument be "time before present" 
#	in getRateThroughTimeMatrix
#


