\name{summary.distfun}
\alias{summary.distfun}
\alias{summary.funxy}
\title{
  Summarizing a Function of Spatial Location
}
\description{
  Prints a useful summary of a function of spatial location.
}
\usage{
  \method{summary}{distfun}(object, \dots)

  \method{summary}{funxy}(object, \dots)
}
\arguments{
  \item{object}{
    An object of class \code{"distfun"} or \code{"funxy"}
    representing a function of spatial coordinates.
  }
  \item{\dots}{
    Arguments passed to \code{\link{as.mask}} controlling the
    pixel resolution used to compute the summary.
  }
}
\details{
  These are the \code{\link[base]{summary}} methods
  for the classes \code{"funxy"} and \code{"distfun"}.
  
  An object of class \code{"funxy"} represents a function
  of spatial location, defined in a particular region of space.
  This includes objects of the special class \code{"distfun"} which represent
  distance functions.

  The \code{summary} method computes a summary of the function values.
  The function is evaluated on a grid of locations using
  \code{\link{as.im}} and numerical values at these locations
  are summarised using \code{\link{summary.im}}. The pixel resolution
  for the grid of locations is determined by the arguments
  \code{\dots} which are passed to \code{\link{as.mask}}.
}
\value{
  For \code{summary.funxy} the result is
  an object of class \code{"summary.funxy"}.
  For \code{summary.distfun} the result is
  an object of class \code{"summary.distfun"}.
  There are \code{print} methods for these classes.
}
\author{
  \spatstatAuthors.
}
\seealso{
  \code{\link{distfun}}, \code{\link{funxy}}
}
\examples{
   f <- function(x,y) { x^2 + y^2 - 1} 
   g <- funxy(f, square(2))
   summary(g)

   summary(distfun(cells))
   summary(distfun(cells), dimyx=256)
}
\keyword{spatial}
\keyword{methods}
