\name{computeBayesFactors}
\alias{computeBayesFactors}
\title{
computeBayesFactors
}
\description{
Computes pairwise Bayes factors for a set of macroevolutionary models sampled using \code{BAMM}, using MCMC simulation output.
}
\usage{
computeBayesFactors(postdata, priordata, burnin = 0.1, modelset = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{postdata}{
Filename for the MCMC output file from a \code{BAMM} run. Alternatively, a dataframe containing this information.
}
  \item{priordata}{
Filename for the \code{BAMM}-generated output file containing a simulation of shifts under the prior. Alternatively, a dataframe containing this information.
}
  \item{burnin}{
What fraction of samples to discard from postdata as burnin?
}
  \item{modelset}{
Integer set of models for which you wish to compute Bayes factors. An integer value of k implies a model with k rate shifts.
}
 \item{...}{
 Additional arguments to computeBayesFactors
 }
}
\details{
This function returns a matrix of pairwise Bayes factors, where the Bayes factor is the ratio of marginal likelihoods between two models M_{i} and M_{j}. Numerator models are given as rows, and denominator models as columns. Row names and column names give the number of shifts in the corresponding model. Suppose you have an output matrix with row and column names 0:3 (0, 1, 2, 3). Model 0 is a model with just a single process (starting at the root), and no among-lineage rate heterogeneity. 

If \code{computeBayesFactors} gives a matrix \code{mm}, and \code{mm[2,1]} is 10.0, this implies Bayes factor evidence of 10 in favor of the 2nd row model (a model with 1 process; e.g., \code{rownames(mm)[2]} over the first column model (a model with a single process).

It can be very difficult to approximate either the posterior or prior probability of models that are very rarely sampled. Suppose you have a simulation of a prior distribution that samples the models (M0, M1, M2, M3). Now, suppose your posterior distribution samples the models (M4, M5, M6, M7). There is no overlap between these sets of models. In other words, models sampled during the simulation of the prior are \emph{so improbable} that they were not sampled during simulation of the posterior. Hence, we can't even estimate their posterior probabilities with any degree of accuracy. Likewise, the models sampled during simulation of the posterior are very unlikely under the prior, and we can't estimate their prior probabilities with any degree of accuracy. The function \code{computeBayesFactors} will automatically check and make sure that it is only computing Bayes factors for models where a non-zero probability (posterior or prior) can be estimated. 

This is likely to change in the future with more robust methods for estimating these quantities in the tails of the distributions. However, keep in mind that if your posterior is very far from your prior, such that these models don't overlap, you actually have strong evidence in favor of rate shifts (if the mean of the posterior distribution is greater than the prior mean) - this implies that the models expected under the prior alone have effectively zero posterior probability.

}
\value{
A matrix of pairwise Bayes factors between models.
}

\author{

Dan Rabosky

}

\examples{
data(mcmc.whales)
data(prior.whales)
computeBayesFactors(mcmc.whales, prior.whales, burnin=0.1)
}

\keyword{ models }
