\name{corBAMM}
\alias{corBAMM}
\alias{corMatrix.corBAMM}
\alias{coef.corBAMM}
\title{
BAMM correlation structure
}
\description{
Implements a mixture of a phylogenetic correlation matrix and a BAMM correlation matrix
for Generalized Least Squares modeling.
}
\usage{
corBAMM(value = numeric(0), ephy, form = ~1, fixed = FALSE, type = c("lambda","exponential", "linear"))
\method{corMatrix}{corBAMM}(object, covariate = getCovariate(object), corr = TRUE, ...)
\method{coef}{corBAMM}(object, unconstrained = TRUE, ...)
}
\arguments{
  \item{value}{
a numeric vector of length 1 or length 2. The first element is the mixture coefficent.
The second element is Pagel's lambda or a shape parameter depending on the value for 
\code{type}.
}
  \item{ephy}{
an object of class "bammdata".
}
  \item{form}{
ignored.  	
}
  \item{fixed}{
a logical indicating whether the parameters should be estimated or treated as
fixed values.  	
}
  \item{object}{
an initialized object of class "corBAMM".  	
}
  \item{covariate}{
ignored. 	
}
  \item{corr}{
a logical indicating whether to return the correlation or the covariance matrix. 
For corBAMM only the correlation matrix is relevant.   	
}
  \item{type}{
a character string indicating the type of mixture.  	
}
  \item{\dots}{
additional arguments to other methods.  	
} 
}
\details{
The BAMM correlation structure is a linear mixture of two correlation matrices, with the contribution of each 
correlation matrix determined by the mixture coefficient "psi" 

Em := psi*Vp + (1-psi)*Vb, 
where Vp is the phylogenetic correlation matrix and Vb is the BAMM correlation matrix. 

The BAMM correlation matrix, Vb, is determined from the posterior distribution and represents the 
probability that any two taxa share a common macroevolutionary rate regime. It is calculated with the
function \code{getCohortMatrix}.

In addition to estimating psi, the user has the option to transform Vp (through Pagel's lambda)
or to transform Vb (using both linear and exponential transformations).

Lambda transform:
    Vp = lambda*Vp (off-diagonal elements only)
Exponential transform:
    Vb = exp(-alpha*(1-Vb))
Linear transform:
    Vb = 1-beta*(1-Vb)

If a transformation is performed the trasformation coefficients will be estimated in addition to
the mixture coefficient.

If no value is specified only psi is estimated and neither matrix is transformed. Specifying a single
value is equivalent to specifying no value, it only changes the initial starting value for psi. If a second value is specified this will be the initial value of the transformation coefficent, which may lambda, alpha, or beta depending
on the type of transformation. Note that if fixed = FALSE all coefficients will remain at their starting
values.
}
\value{
an object of class "corBAMM", the coeffients of an object of this class, or the correlation matrix of 
an initialized object of this class.  In most cases the user will only need to call \code{corBAMM}.
}
\author{
Mike Grundler
}
\note{
The exponential and linear transformations of the BAMM correlation matrix are experimental and
may cause difficulties for inference.
}

\seealso{
\code{\link{getCohortMatrix}}, \code{\link{vcv.phylo}}, \code{\link{corPagel}}
}
\examples{
 
}
\keyword{ model }
