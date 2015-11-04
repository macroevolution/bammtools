\name{plotPrior}
\alias{plotPrior}
\title{
Plot the prior and posterior distribution of shifts
}
\description{
Generates a barplot of the prior and posterior distributions of the number of shifts.
}
\usage{
plotPrior(mcmc, expectedNumberOfShifts = 1, burnin = 0.15, 
	priorCol = 'light blue', postCol = 'red', legendPos = 'topright', ...)
}
\arguments{
	\item{mcmc}{
		a dataframe of the mcmc_out file from a BAMM run, or the filename
	}
	\item{expectedNumberOfShifts}{
		Expected number of shifts under the prior
	}
	\item{burnin}{
		the fraction of samples to discard as burn-in
	}
	\item{priorCol}{
		color for the prior distribution
	}
	\item{postCol}{
		color for the posterior distribution
	}
	\item{legendPos}{
		placement of the legend, see \code{\link{legend}}
	}
	\item{...}{
		additional parameters that are passed to \code{\link{barplot}}
	}
}


\value{
Invisibly returns a matrix with the probability of each shift number
under the prior and the posterior.


}
\author{
	Pascal Title
}
\examples{
data(mcmc.whales)

plotPrior(mcmc.whales, expectedNumberOfShifts = 1, burnin = 0.15)
}
 

