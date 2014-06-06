\name{distinctShiftConfigurations}
\alias{distinctShiftConfigurations}

\title{
Identify distinct rate shift configurations
}
\description{
Identify topologically distinct rate shift configurations that were sampled with BAMM, and assign each sample in the posterior to one of the distinct shift configurations.	
}
\usage{
distinctShiftConfigurations(ephy, prior, BFcriterion, ...)
}
\arguments{
	\item{ephy}{
		an object of class \code{bammdata}.
	}
	\item{prior}{
		an object of class \code{branchprior}
	}
	
	\item{BFcriterion}{
		Numeric value for the Bayes factor threshold used to identify branches with elevated shift probabilities relative to prior
	}
	\item{...}{
		Other arguments to distinctShiftConfigurations (possibly deprecated args)
	}
}

\value{
an object of class \code{bammshifts}. This is a list with the following components
\item{marg.probs}{
A list of the marginal probability of a shift occurring at each node of the phylogeny for each
distinct rate configuration.	
}
\item{bayesfactors}{
Bayes factor evidence for a rate shift on each branch
}
\item{shifts}{
A list of the set of shift nodes for each distinct rate configuration.	
}
\item{samplesets}{
A list of sample indices that reduce to each of the unique shift sets.	
}
\item{frequency}{
A vector of frequencies of each distinct shift configuration	.
}
\item{coreshifts}{
A vector of node numbers corresponding to the core shifts. All of these nodes have a Bayes factor of at least \code{BFcriterion} supporting a rate shift.
}
\item{BFcriterion}{
A single numeric value giving the Bayes factor threshold used during enumeration of distinct shift configurations
}	
Results are sorted by frequency:
$frequency[1] gives the most common shift configuration sampled.
$shifts[[1]] gives the corresponding node indices for that configuration.
$samplesets[[1]] gives the indices of samples with this configuration.
}

\author{
	Dan Rabosky
}
\seealso{
	\code{\link{plot.bammshifts}}, \code{\link{credibleShiftSet}}
}
\examples{
data(whales, events.whales, prior.whales)

# Get the prior distribution on the number of shifts per branch:
bp <- getBranchShiftPriors(whales, prior.whales)


ed <- getEventData(whales, events.whales, burnin=0.25)

sc <- distinctShiftConfigurations(ed, bp, 10)

plot.bammshifts(sc, ed, rank=1)

}



