\name{layered}
\alias{layered}
\title{
  Create List of Plotting Layers
}
\description{
  Given several objects which are capable of being plotted,
  create a list containing these objects as if they were
  successive layers of a plot. The list can then be plotted
  in different ways.
}
\usage{
layered(..., plotargs = NULL, LayerList=NULL)
}
\arguments{
  \item{\dots}{
    Objects which can be plotted by \code{plot}.
  }
  \item{plotargs}{
    Default values of the plotting arguments for each of the objects.
    A list of lists of arguments of the form \code{name=value}.
  }
  \item{LayerList}{
    A list of objects.
    Incompatible with \code{\dots}.
  }
}
\details{
  Layering is a simple mechanism for controlling
  a high-level plot that is composed of
  several successive plots, for example, a background and a foreground
  plot. The layering mechanism makes it easier to issue the plot command,
  to switch on or off the plotting of each individual layer,
  to control the plotting arguments that are passed to each layer,
  and to zoom in.
  
  Each individual layer in the plot should be saved as an object
  that can be plotted using \code{plot}. It will typically belong to
  some class, which has a method for the generic function \code{plot}.

  The command \code{layered} simply saves the objects \code{\dots}
  as a list of class \code{"layered"}. This list can then be plotted by
  the method \code{\link{plot.layered}}. Thus, you only need to
  type a single \code{plot} command to produce the multi-layered plot.
  Individual layers of the plot can be switched on or off, or
  manipulated, using arguments to \code{\link{plot.layered}}.

  The argument \code{plotargs} contains default values of the
  plotting arguments for each layer. It should be a list, with one
  entry for each object in \code{\dots}. Each entry of  \code{plotargs}
  should be a list of arguments in the form \code{name=value}, which are
  recognised by the \code{plot} method for the relevant layer.

  The \code{plotargs} can also include an argument named \code{.plot}
  specifying (the name of) a function to perform the plotting
  instead of the generic \code{plot}.

  The length of \code{plotargs} should either be equal to the
  number of layers, or equal to 1. In the latter case it will be
  replicated to the appropriate length.
}
\value{
  A list, belonging to the class \code{"layered"}.
  There are methods for \code{plot}, \code{"["},
  \code{"shift"}, \code{"affine"}, \code{"rotate"} and \code{"rescale"}.
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\seealso{
  \code{\link{plot.layered}},
  \code{\link{methods.layered}},
  \code{\link{as.layered}},
  \code{\link{[.layered}},
  \code{\link{layerplotargs}}.
}
\examples{
   D <- distmap(cells)
   L <- layered(D, cells)
   L
   L <- layered(D, cells,
    plotargs=list(list(ribbon=FALSE), list(pch=16)))
   plot(L)

   layerplotargs(L)[[1]] <- list(.plot="contour")
   plot(L)
}
\keyword{spatial}
\keyword{hplot}
