\name{pairs.linim}
\alias{pairs.linim}
\title{
  Scatterplot Matrix for Pixel Images on a Linear Network
}
\description{
  Produces a scatterplot matrix of the pixel values
  in two or more pixel images on a linear network.
}
\usage{
\method{pairs}{linim}(..., plot=TRUE, eps=NULL)
}
\arguments{
  \item{\dots}{
    Any number of arguments, each of which is either
    a pixel image on a linear network (object of class \code{"linim"}),
    a pixel image (object of class \code{"im"}),
    or a named argument to be passed to \code{\link{pairs.default}}.
  }
  \item{plot}{
    Logical. If \code{TRUE}, the scatterplot matrix is plotted.
  }
  \item{eps}{
    Optional. Spacing between sample points on the network.
    A positive number.
  }
}
\details{
  This is a method for the generic function \code{\link{pairs}}
  for the class of pixel images on a linear network.
  
  It produces a square array of plot panels, in which each panel shows
  a scatterplot of the pixel values of one image against the
  corresponding pixel values of another image.

  At least two of the arguments \code{\dots} should be a pixel image
  on a linear network (object of class \code{"linim"}).
  They should be defined on the \bold{same} linear network,
  but may have different pixel resolutions.

  First the pixel values of each image are extracted at a
  set of sample points equally-spaced across the network.
  Then \code{\link{pairs.default}} is called to 
  plot the scatterplot matrix.

  Any arguments in \code{\dots} which are not pixel images will be
  passed to \code{\link{pairs.default}} to control the plot.
}
\value{
  Invisible. A \code{data.frame} containing the
  corresponding pixel values for each image.
  The return value also belongs to the class \code{plotpairsim} which has
  a plot method, so that it can be re-plotted.
}
\seealso{
  \code{\link{pairs.default}},
  \code{\link{pairs.im}}
}
\examples{
  fit <- lppm(chicago ~ marks * (x+y))
  lam <- predict(fit)
  do.call(pairs, lam)
}
\author{
  \spatstatAuthors
}
\keyword{spatial}
\keyword{hplot}

