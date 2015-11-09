\name{addBAMMlegend}
\alias{addBAMMlegend}
\title{
Add a color legend to a phylo-rate plot
}
\description{
Add a legend to a phylorate plot, with greater manual control.
}
\usage{
addBAMMlegend(x, direction, side, location = 'topleft', nTicks = 2, shortFrac = 0.02, 
longFrac = 0.3, axisOffset = 0.002, cex.axis = 0.8, labelDist = 0.7, ...)
}
\arguments{
	\item{x}{
		a plot.bammdata object.
	}
	\item{direction}{
		direction of color ramp. If omitted, then direction is automatically inferred,
		otherwise can be specified as \code{horizontal} or \code{vertical}.
	}
	\item{side}{
		side for tick marks, see \code{\link{axis}} documentation. Automatically inferred if omitted.
	}
	\item{location}{
		either a location name (see \code{Details}), or coordinates for the corners of the bar legend \code{c(xmin, xmax, ymin, ymax)}.	
	}
	\item{nTicks}{
		number of tick marks, besides of min and max.
	}
	\item{shortFrac}{
		Percent of the plot width range that will be used as the short dimention of the legend. Only applies to preset location options. 
	}
	\item{longFrac}{
		Percent of the plot width range that will be used as the long dimention of the legend. Only applies to preset location options.	
	}
	\item{axisOffset}{
		distance from color bar for labels, as a percent of the plot width.
	}
	\item{cex.axis}{
		size of axis labels
	}
	\item{labelDist}{
		distance from axis to axis labels (passed to \code{mgp})
	}
	\item{...}{
		additional parameters to be passed to \code{\link{axis}}.
	}
}

\details{
	A number of predefined locations exist in this function to make it easy to add a legend to a phylorate plot. 
	Preset \code{locations} are: \code{topleft}, \code{topright}, \code{bottomleft}, \code{bottomright}, 
	\code{left}, \code{right}, \code{top} and \code{bottom}. 
	If more fine-tuned control is desired, then a numeric vector of length 4 can be supplied to \code{location}, 
	specifying the min x, max x, min y and max y values for the legend. See examples.
}

\value{
	Invisibly returns a list with the following components.

	\item{coords}{
		2-column matrix of xy coordinates for each color bin in the legend. 
	}
	\item{width}{
		Coordinates for the short dimension of the legend. 
	}
	\item{pal}{
		the color ramp
	}
	\item{tickLocs}{
		the tick mark locations in plotting units
	}
	\item{labels}{
		the rate values associated with those tick locations. 
	}
}
\author{
	Pascal Title
}
\seealso{
	Requires an object created with \code{\link{plot.bammdata}}. 
}
\examples{
data(whales, events.whales)
ephy <- getEventData(whales, events.whales, burnin = 0.25, nsamples = 300)

# plot phylorate with extra margin space
x <- plot.bammdata(ephy, lwd = 2, mar = c(5,4,4,4)) 
# presets
addBAMMlegend(x, location = 'topleft')
addBAMMlegend(x, location = 'bottom')
addBAMMlegend(x, location = 'right')

#fine-tune placement
x <- plot.bammdata(ephy, lwd = 2, mar = c(5,4,4,4)) 
axis(1); axis(2)
addBAMMlegend(x, location = c(-1, -0.5, 40, 80), nTicks = 4)
addBAMMlegend(x, location = c(5, 20, 60, 61), nTicks = 4, side = 3, cex.axis = 0.7)

#addBAMMlegend also automatically detects the use of color.interval
data(primates, events.primates)
ephy <- getEventData(primates, events.primates, burnin=0.25, nsamples = 300, type = 'trait')

x <- plot.bammdata(ephy, breaksmethod = 'linear', color.interval = c(NA, 0.12), lwd = 2)
addBAMMlegend(x, location = c(0, 30, 200, 205), nTicks = 1, side = 3)
}
 

