data(whales, events.whales, prior.whales)

# Get the prior distribution on the number of shifts per branch:
bp <- getBranchShiftPriors(whales, prior.whales)


ed <- getEventData(whales, events.whales, burnin=0.25)

sc <- distinctShiftConfigurations(ed, bp, 10)

plot.bammshifts(sc, ed, rank=1)