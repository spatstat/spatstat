\name{contour.listof}
\alias{contour.listof}
\alias{image.listof}
\title{Plot a List of Things}
\description{
  Plots a list of things, either as an array of contour plots,
  or as an array of images.
}
\usage{
  \method{contour}{listof}(x, \dots)
  \method{image}{listof}(x, \dots, equal.ribbon=FALSE)
}
\arguments{
  \item{x}{
    An object of the class \code{"listof"}.
    Essentially a list of objects.
  }
  \item{\dots}{
    Arguments passed to \code{\link{plot.listof}} to control the plot.
  }
  \item{equal.ribbon}{
    Logical. If \code{TRUE}, the colour maps of all the images are
    the same. If \code{FALSE}, the colour map of each
    image is adjusted to the range of values of that image.
  }
}
\value{
  Null.
}
\details{
  These are methods for the generic plot commands
  \code{contour} and \code{image}, for the class \code{"listof"}.

  The commands will display each object in the list \code{x},
  either as a contour plot
  (\code{contour.listof}) or a pixel image plot (\code{image.listof}),
  with the displays laid out in a grid. 
  See \code{\link{plot.listof}} for more information.
}
\seealso{
  \code{\link{plot.listof}}
}
\examples{
# Multitype point pattern
 contour(D <- density(split(amacrine)))
 image(D, equal.ribbon=TRUE, main="")
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{hplot}
