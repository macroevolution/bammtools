\name{assignColorBreaks}
\alias{assignColorBreaks}
\title{
	Map macroevolutionary rates to colors
}
\description{
	Maps macroevolutionary rates to a set of \code{NCOLORS}.
}
\usage{
	assignColorBreaks(rates, NCOLORS, spex, log = FALSE, method = c("linear","quantile","jenks"), 
		color.interval = NULL, JenksSubset = NULL)
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
	\item{color.interval}{
		Numeric vector of length 2. Min and max values for the linear mapping of rates to colors. 
		If \code{NULL}, then min and max are inferred from the data. Only relevant 
		when \code{method = "linear"}. See Details. 
	}
	\item{JenksSubset}{
		Number of regularly spaced samples to subset from \code{rates}. Only relevant when \code{method = "jenks"}. 
		See Details. 
	}
}
\details{
	If \code{method = "quantile"} macroevolutionary rates are binned into \code{NCOLORS+1} percentiles 
	and rates in each bin are mapped to a color determined by the \code{pal} argument in \code{plot.bammdata}. 
	Alternatively, if \code{method = "linear"} macroevolutionary rates are binned into \code{NCOLORS+1} equal 
	length intervals between the minimum and maximum. If color.interval is defined, then those min and max values 
	override the automatic detection of min and max. This might be useful if some small number of lineages have very
	high or very low rates, such that the linear map of colors is being skewed towards these extremes, resulting in
	other rate variation being drowned out. 
	If \code{method = "jenks"}, macroevolutionary rates are binned into \code{NCOLORS+1} categories, according to 
	the Jenks natural breaks classification method. This method is borrowed from the field of cartography, and 
	seeks to minimize the variance within categories, while maximizing the variance between categories. 

	The Jenks natural breaks method was ported to C from code found in the classInt R package. 
}
\value{
	A numeric vector of rate percentiles/intervals. 
}
\author{
	Mike Grundler, Pascal Title
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