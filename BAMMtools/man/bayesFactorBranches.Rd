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
 \item{priordata}{
an object of class \code{branchprior}
 }
}
\details{
This function returns a copy of a phylogenetic tree where each branch length is equal to the Bayes factor evidence in favor of a rate shift on a particular branch. Here, we imagine that there are two models: one with a rate shift on a given branch, and one without. Treating the "no shift" model as something akin to a null model, the Bayes factor can be computed as the ratio of posterior to prior odds of the model with a shift on a given branch, divided by the ratio of posterior to prior odds of the model without a shift. 

By computing branch-specific marginal shift probabilities (using \code{\link{marginalShiftProbsTree}}), you can immediately extract the posterior probabilities of the relevant models. The branch length returned by this function is the posterior probability of a rate shift (\emph{prob_shift}), and the posterior probability of the model with no shift is simply \emph{1 - prob_shift}. 

Likewise, you can compute the prior odds in much the same way. The function \code{\link{getBranchShiftPriors}} computes the \emph{prior} probability of a rate shift on each branch. As above, the prior probability of the \emph{no shift} model is just \emph{1 - prob_noshift}. The posterior-to-prior odds ratio for the shift model is \emph{(posterior prob of shift) / (prior prob of shift)}. Simply multiplying this term by \emph{(prior prob of no shift) / (posterior prob of no shift)} gives us the Bayes factor. 

In general, Bayes factors greater than 10 imply evidence in support of a rate shift, and Bayes factors greater than 100 imply very strong evidence for a shift on a particular branch. 
 
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
plot.phylo(whales, show.tip.label=FALSE)
mtext(side=3, text="original tree", cex=1)

plot.phylo(margprobs, show.tip.label=FALSE)
mtext(side=3, text="Marginal probs", cex=1)

plot.phylo(bftree, show.tip.label=FALSE)
mtext(side=3, text="Bayes factors", cex=1)


}
