\name{bayesFactorBranches}
\alias{bayesFactorBranches}
\title{
Bayes factor evidence for rate shifts on individual branches
}
\description{
Computes Bayes factor evidence for a rate shift on every branch of a phylogenetic tree
}
\usage{
bayesFactorBranches(ephy, priordata)
}
\arguments{
  \item{ephy}{
an object of class \code{bammdata}
}
 \item{prior}{
an object of class \code{branchprior}
 }
}
\details{
This function returns a copy of a phylogenetic tree where each branch length is equal to the Bayes factor evidence in favor of a rate shift on a particular branch. The Bayes factor has a natural interpretation in this context: it is simply the ratio of the posterior probability of a rate shift on a given branch normalized by the prior probability of a rate shift on the same branch. 

A branch length of 1 implies that the observed posterior shift probability is identical to that expected under the prior alone. A branch length of 100 implies that the posterior probability of a shift is 100 times higher than the prior. A value of 0.1 implies that a shift is actually less likely on a particular branch than expected under the prior alone.
}

\value{
A object of class \code{phylo} but where each branch length is equal to the Bayes factor evidence in favor of a rate shift on a particular branch. 
}

\author{
Dan Rabosky
}

\seealso{
\code{\link{getBranchShiftPriors}}, \code{\link{distinctShiftConfigurations}}, \code{\link{credibleShiftSet}}
}
\examples{
data(whales, events.whales, prior.whales)
ed <- getEventData(whales, events.whales, burnin=0.1)
pdata <- getBranchShiftPriors(whales, prior.whales)
bftree <- bayesFactorBranches(ed, pdata)
margprobs <- marginalShiftProbsTree(ed)
plot.new()
par(mfcol=c(1,3))
# Original tree:
plot.phylo(whales, show.tip.label=F)
mtext(side=3, text="original tree", cex=1)

plot.phylo(margprobs, show.tip.label=F)
mtext(side=3, text="Marginal probs", cex=1)

plot.phylo(bftree, show.tip.label=F)
mtext(side=3, text="Bayes factors", cex=1)


}
