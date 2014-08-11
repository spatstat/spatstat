\name{methods.objsurf}
\Rdversion{1.1}
\alias{methods.objsurf} %DoNotExport
\alias{print.objsurf}
\alias{plot.objsurf}
\alias{persp.objsurf}
\alias{image.objsurf}
\alias{contour.objsurf}
\title{
  Methods for Objective Function Surfaces
}
\description{
  Methods for printing and plotting an objective function surface.
}
\usage{
\method{print}{objsurf}(x, ...)
\method{plot}{objsurf}(x, ...)
\method{image}{objsurf}(x, ...)
\method{contour}{objsurf}(x, ...)
\method{persp}{objsurf}(x, ...)
}
\arguments{
  \item{x}{
    Object of class \code{"objsurf"} representing an objective function surface.
  }
  \item{\dots}{
    Additional arguments passed to plot methods.
  }
}
\details{
  These are methods for the generic functions
  \code{\link{print}}, \code{\link{plot}},
  \code{\link{image}}, \code{\link{contour}} and  \code{\link{persp}}
  for the class \code{"objsurf"}.
}
\value{
  For \code{print.objsurf}, \code{plot.objsurf} and \code{image.objsurf}
  the value is \code{NULL}.

  For \code{contour.objsurf} and \code{persp.objsurf}
  the value is described in the help for
  \code{\link{contour.default}} and \code{\link{persp.default}}
  respectively.
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Ege Rubak.
}
\seealso{
  \code{\link{objsurf}}
}
\examples{
  example(objsurf)
  os
  plot(os)
  contour(os, add=TRUE)
  persp(os)
}
\keyword{spatial}
\keyword{hplot}
