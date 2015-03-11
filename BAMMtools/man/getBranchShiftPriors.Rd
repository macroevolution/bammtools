\name{getBranchShiftPriors}
\alias{getBranchShiftPriors}
\title{
Compute prior odds of a rate shift on each branch of a phylogeny from BAMM output
}
\description{
Computes the prior probability of a rate shift event for each branch. These results are important for identifying topological rate shift locations on phylogenies with marginal probabilities that exceed those predicted under the prior alone. 
}
\usage{
getBranchShiftPriors(phy, priordata)
}
\arguments{
  \item{phy}{
an object of class \code{phylo}
}
  \item{priordata}{
Filename for the \code{BAMM}-generated output file containing a simulation of shifts under the prior. Alternatively, a dataframe containing this information. 
}
}
\details{

This function computes the prior odds on the distribution of numbers of rate shift events per branch under the prior. It returns an object of class 'branchprior', which is nothing more than a copy of the original phylogenetic tree but where each branch length has been replaced by the prior probability of a rate shift on each branch.

The significance of this function is that it lets us explicitly determine which branches have shift probabilities that are elevated relative to the prior expectation. The user does not work with the output of this function directly; it is passed to other functions \link[BAMMtools]{credibleShiftSet} and \link[BAMMtools]{distinctShiftConfigurations}.

Note: If you are working with a subset of a bammdata object (created with \code{\link{subtreeBAMM}}), then you must generate a new branchpriors object with a pruned phylogeny that matches the one associated with your subsetted bammdata object. 
 
}
\value{
A class \code{branchprior} with all the components of the original class \code{phylo} object, with the following changes:

\item{edge.length}{
branch lengths now represent the prior probability of a rate shift on each branch
}
 
}
\references{
\code{bamm-project.org}
}
\author{
Dan Rabosky
}

\seealso{
\code{\link{distinctShiftConfigurations}}, \code{\link{plot.bammshifts}}, \code{\link{summary.credibleshiftset}}, \code{\link{plot.credibleshiftset}}, \code{\link{credibleShiftSet}}	
} 
 
\examples{
data(events.whales)
data(whales)
data(prior.whales)
ed <- getEventData(whales, events.whales)

# Get branch-specific cutoffs for identifying branches with 
#	elevated shift probabilities:
bprior <- getBranchShiftPriors(whales, prior.whales)
cset <- credibleShiftSet(ed, bprior)

# Here is the total number of samples in the posterior:
length(ed$eventData)

# And here is the number of distinct shift configurations:
cset$number.distinct

# here is the summary statistics:
summary(cset)

# Accessing the raw frequency vector for the credible set:
cset$frequency

#The cumulative frequencies:
cset$cumulative

# The first element is the shift configuration with the maximum 
#  a posteriori probability. We can identify all the samples from 
#  posterior that show this shift configuration:

cset$indices[[1]]
 
 
 
}
 
\keyword{ models }
