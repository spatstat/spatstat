\name{pairs.im}
\alias{pairs.im}
\title{
  Scatterplot Matrix for Pixel Images
}
\description{
  Produces a scatterplot matrix of the pixel values
  in two or more pixel images.
}
\usage{
\method{pairs}{im}(..., plot=TRUE)
}
\arguments{
  \item{\dots}{
    Any number of arguments, each of which is either
    a pixel image (object of class \code{"im"})
    or a named argument to be passed to \code{\link{pairs.default}}.
  }
  \item{plot}{
    Logical. If \code{TRUE}, the scatterplot matrix is plotted.
  }
}
\details{
  This is a method for the generic function \code{\link{pairs}}
  for the class of pixel images.
  
  It produces a square array of plot panels, in which each panel shows
  a scatterplot of the pixel values of one image against the
  corresponding pixel values of another image.

  At least two of the arguments \code{\dots}
  should be pixel images (objects of class \code{"im"}).
  Their spatial domains must overlap, but need not have the same
  pixel dimensions.

  First the pixel image domains are intersected, and converted to a
  common pixel resolution. Then the corresponding pixel values of each image are
  extracted. Then \code{\link{pairs.default}} is called to 
  plot the scatterplot matrix.

  Any arguments in \code{\dots} which are not pixel images will be
  passed to \code{\link{pairs.default}} to control the plot.
}
\note{
  To control the appearance of the individual scatterplot panels,
  see \code{\link{pairs.default}}, \code{\link{points}}
  or \code{\link{par}}.
  To control the plotting symbol for the points in the scatterplot,
  use the arguments \code{pch}, \code{col}, \code{bg} as described
  under \code{\link{points}} (because 
  the default panel plotter is the function \code{\link{points}}).
  To suppress the tick marks on the plot axes,
  type \code{par(xaxt="n", yaxt="n")} before calling \code{pairs}.
}
\value{
  Invisible. A \code{data.frame} containing the
  corresponding pixel values for each image.
  The return value also belongs to the class \code{plotpairsim} which has
  a plot method, so that it can be re-plotted.
}
\seealso{
  \code{\link{pairs}},
  \code{\link{pairs.default}},
  \code{\link{plot.im}},
  \code{\link{im}},
  \code{\link{par}}
}
\examples{
  X <- density(rpoispp(30))
  Y <- density(rpoispp(40))
  Z <- density(rpoispp(30))
  pairs(X,Y,Z)
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{hplot}

