\name{Extract.owin}
\alias{[.owin}
\title{Extract Subset of Window}
\description{
  Extract a subset of a window.
}
\usage{
  \method{[}{owin}(x, i, \dots)
}
\arguments{
  \item{x}{
    A spatial window (object of class \code{"owin"}).
  }
  \item{i}{
    Object defining the subregion.
    Either a spatial window, or a
    pixel image with logical values.
  }
  \item{\dots}{Ignored.}
}
\value{
  Another spatial window (object of class \code{"owin"}).
}
\details{
  This function computes the intersection
  between the window \code{x} and the domain specified by \code{i},
  using \code{\link{intersect.owin}}.
  
  This function is a method for the subset operator \code{"["} for
  spatial windows (objects of class \code{"owin"}). It is provided
  mainly for completeness.

  The index \code{i} may be either a window, or a pixel image with
  logical values (the \code{TRUE} values of the
  image specify the spatial domain).
}
\seealso{
  \code{\link{intersect.owin}}
}
\examples{
 W <- owin(c(2.5, 3.2), c(1.4, 2.9))
 plot(letterR)
 plot(letterR[W], add=TRUE, col="red")
}
\author{
  \adrian
  and
  \rolf
}
\keyword{spatial}
\keyword{manip}
