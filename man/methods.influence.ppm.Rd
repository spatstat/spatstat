\name{methods.influence.ppm}
\alias{methods.influence.ppm} %DoNotExport
\alias{as.ppp.influence.ppm}
\alias{as.owin.influence.ppm}
\alias{domain.influence.ppm}
\alias{integral.influence.ppm}
\alias{Smooth.influence.ppm}
\alias{Window.influence.ppm}
\title{Methods for Influence Objects}
\description{
  Methods for the class \code{"influence.ppm"}.
}
\usage{
  \method{as.ppp}{influence.ppm}(X, \dots)

 \method{as.owin}{influence.ppm}(W, \dots, fatal=TRUE)

 \method{domain}{influence.ppm}(X, \dots)

 \method{Smooth}{influence.ppm}(X, \dots)

 \method{Window}{influence.ppm}(X, \dots)

 \method{integral}{influence.ppm}(f, domain, \dots)
}
\arguments{
  \item{X,W,f}{An object of class \code{"influence.ppm"}.}
  \item{domain}{
    Optional. Domain of integration: a window (class \code{"owin"})
    or a tessellation (class \code{"tess"}).
  }
  \item{\dots}{
    Additional arguments.
    See Details.
  }
  \item{fatal}{
    Logical value indicating what to do if the data cannot be
    converted to a window. If \code{fatal=TRUE} (the default)
    an error occurs. If \code{fatal=FALSE} a value of \code{NULL}
    is returned.
  }
}
\value{
  A window (object of class \code{"owin"}) for
  \code{as.owin}, \code{domain} and \code{Window}.
  A point pattern (object of class \code{"ppp"}) for
  \code{as.ppp}. A numeric value or numeric vector for \code{integral}.
  A pixel image, or list of pixel images, for \code{Smooth}.
}
\details{
  These functions are methods for the class \code{"influence.ppm"}.
  An object of this class
  represents the influence measure of a fitted point process model
  (see \code{\link{influence.ppm}}).

  For \code{as.ppp}, \code{domain}, \code{integral} and \code{Window},
  additional arguments (\code{\dots}) are ignored.
  For \code{as.owin} and \code{Smooth}, 
  additional arguments are passed to the method for class \code{"ppp"}.
}
\seealso{
  \code{\link{influence.ppm}}, \code{\link{plot.influence.ppm}},
  \code{\link{[.influence.ppm}}
  
}
\examples{
 fit <- ppm(cells ~ x)
 a <- influence(fit)
 Window(a)
}
\author{
  \spatstatAuthors
}
\keyword{spatial}
\keyword{manip}
