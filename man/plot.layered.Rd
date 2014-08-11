\name{plot.layered}
\alias{plot.layered}
\title{
Layered Plot
}
\description{
  Generates a layered plot.
  The plot method for objects of class \code{"layered"}.
}
\usage{
\method{plot}{layered}(x, ..., which = NULL, plotargs = NULL)
}
\arguments{
  \item{x}{
    An object of class \code{"layered"}
    created by the function \code{\link{layered}}.
  }
  \item{\dots}{
    Arguments to be passed to the \code{plot} method
    for \emph{every} layer.
  }
  \item{which}{
    Subset index specifying which layers should be plotted. 
  }
  \item{plotargs}{
    Arguments to be passed to the  \code{plot} methods
    for individual layers. A list of lists of arguments of the
    form \code{name=value}.
  }
}
\details{
  Layering is a simple mechanism for controlling
  a high-level plot that is composed of
  several successive plots, for example, a background and a foreground
  plot. The layering mechanism makes it easier to plot,
  to switch on or off the plotting of each individual layer,
  and to control the plotting arguments that are
  passed to each layer.
  
  The layers of data to be plotted should first be converted
  into a single object of class \code{"layered"} using the
  function \code{\link{layered}}. Then the layers can be plotted using
  the method \code{plot.layered}.

  The subset operator \code{\link{[.layered}} can be used
  to zoom in on a subregion.
}
\value{
  Null.
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\seealso{
  \code{\link{layered}},
  \code{\link{layerplotargs}},
  \code{\link{[.layered}},
  \code{\link{plot}}.
}
\examples{
   data(cells)
   D <- distmap(cells)
   L <- layered(D, cells)
   plot(L)
   plot(L, which = 2)
   plot(L, plotargs=list(list(ribbon=FALSE), list(pch=3, cols="white")))
   # plot a subregion
   plot(L[, square(0.5)])
}
\keyword{spatial}
\keyword{hplot}
