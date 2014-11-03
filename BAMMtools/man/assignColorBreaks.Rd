\name{assignColorBreaks}
\alias{assignColorBreaks}
\title{
	Map macroevolutionary rates to colors
}
\description{
	Maps macroevolutionary rates to a set of \code{NCOLORS}.
}
\usage{
	assignColorBreaks(rates, NCOLORS, spex, log = FALSE, method = c("linear", "quantile"))
}
\arguments{
	\item{rates}{
		a numeric vector of phenotypic rates or a list of numeric vectors of speciation and extinction rates.
	}
	\item{NCOLORS}{
		an integer number of colors to use for the mapping. Larger numbers do not 
		necessarily result in smoother looking color ramps. The default is 64 and 
		is probably sufficient for most purposes.
	}
	\item{spex}{
		a character string. "s" means that speciation rates are used to make the
		map, "e" means that extinction rates are used. Any other string, e.g. "se",
		means that diversification rates are used. Ignored for BAMM trait data.
	}
	\item{log}{
		should the natural logarithm of rates be used for the color map.
	}
	\item{method}{
		determines how the color breaks are created. See details.
	}
}
\details{
	If \code{method = "quantile"} macroevolutionary rates are binned into \code{NCOLORS+1} percentiles and rates in 	each bin are mapped to a color determined by the \code{pal} argument in \code{plot.bammdata}. Alternatively,
	if \code{method = "linear"} macroevolutionary rates are binned into \code{NCOLORS+1} equal length intervals		between the minimum and maximum.
}
\value{
	A numeric vector of rate percentiles/intervals. 
}
\author{
	Mike Grundler
}
\examples{
	data(whales, events.whales)
	ed <- getEventData(whales,events.whales,0.2)
	
	ed <- dtRates(ed, 0.01)
	colors <- assignColorBreaks(ed$dtrates$rates, spex="s") #speciation rates
	#colors <- assignColorBreaks(ed$dtrates$rates[[1]]) 
	#this also works for speciation rates
	
	plot(ed, colorbreaks = colors, spex="s")
	colors <- assignColorBreaks(ed$dtrates$rates, spex="se") 
	#diversification rates
	
	#colors <- assignColorBreaks(ed$dtrates$rates[[1]] - ed$dtrates$rates[[2]]) 
	#this also works for diversification rates
	
	plot(ed, colorbreaks = colors, spex="se")
}
\seealso{
	\code{\link{plot.bammdata}}
}
\keyword{ graphics }