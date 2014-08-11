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
    Arguments passed to \code{\link{plot.listof}} to control the plot,
    and passed to the appropriate plot command for each element of the list.
  }
  \item{equal.ribbon}{
    Logical. If \code{TRUE}, the colour maps of all the images are
    the same. If \code{FALSE}, the colour map of each
    image is adjusted to the range of values of that image.
    Applies only when all the elements of \code{x} are
    pixel images (objects of class \code{"im"}).
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

  For the image plot, if \code{equal.ribbon=FALSE},
  the images are rendered using different colour maps,
  which are displayed as colour ribbons beside each image.
  If \code{equal.ribbon=TRUE}, the images are rendered using the
  same colour map, and a single colour ribbon will be displayed at the right
  side of the array. 
  The colour maps and the placement of the colour ribbons are
  controlled by arguments \code{\dots} passed to \code{\link{plot.im}}.
}
\seealso{
  \code{\link{plot.listof}},
  \code{\link{plot.im}}
}
\examples{
# Multitype point pattern
 contour(D <- density(split(amacrine)))
 image(D, equal.ribbon=TRUE, main="", col.ticks="red", col.axis="red")
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{hplot}
