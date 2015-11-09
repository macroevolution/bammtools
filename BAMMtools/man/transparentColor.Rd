\name{transparentColor}
\alias{transparentColor}

\title{
Define colors with transparency
}
\description{
Converts a named color and opacity and returns the proper RGB code.
}
\usage{
transparentColor(namedColor, alpha = 0.8)
}

\arguments{
  \item{namedColor}{
a color name
}
  \item{alpha}{
a transparency value between 0 and 1, where 0 is fully transparent
}
}
\details{
This function is used internally by \code{\link{plotRateThroughTime}}. 
}
\value{
Returns the transparent color in RGB format. 
}
 
\author{
Pascal Title
}
 
\keyword{ manip }
