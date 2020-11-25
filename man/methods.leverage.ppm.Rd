\name{methods.leverage.ppm}
\alias{methods.leverage.ppm} %DoNotExport
\alias{as.im.leverage.ppm}
\alias{as.owin.leverage.ppm}
\alias{domain.leverage.ppm}
\alias{integral.leverage.ppm}
\alias{mean.leverage.ppm}
\alias{Smooth.leverage.ppm}
\alias{Window.leverage.ppm}
\title{Methods for Leverage Objects}
\description{
  Methods for the class \code{"leverage.ppm"}.
}
\usage{
  \method{as.im}{leverage.ppm}(X, \dots, what=c("smooth", "nearest"))

 \method{as.owin}{leverage.ppm}(W, \dots, fatal=TRUE)

 \method{domain}{leverage.ppm}(X, \dots)

 \method{integral}{leverage.ppm}(f, domain, \dots)

 \method{mean}{leverage.ppm}(x, \dots)

 \method{Smooth}{leverage.ppm}(X, \dots)

 \method{Window}{leverage.ppm}(X, \dots)

}
\arguments{
  \item{X,x,W,f}{An object of class \code{"leverage.ppm"}.}
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
 \item{what}{
   Character string (partially matched) specifying which image data
   should be extracted. See \code{\link{plot.leverage.ppm}} for
   explanation.
 }
}
\value{
  A window (object of class \code{"owin"}) for
  \code{as.owin}, \code{domain} and \code{Window}.
  A numeric value or numeric vector for \code{integral}.
  A pixel image, or list of pixel images, for \code{as.im} and \code{Smooth}.
}
\details{
  These functions are methods for the class \code{"leverage.ppm"}.
  An object of this class
  represents the leverage measure of a fitted point process model
  (see \code{\link{leverage.ppm}}).

  For \code{as.im}, \code{domain} and \code{Window},
  additional arguments (\code{\dots}) are ignored.
  For \code{as.owin}, \code{integral}, \code{mean} and \code{Smooth}, 
  additional arguments are passed to the method for class \code{"im"}.
}
\seealso{
  \code{\link{leverage.ppm}}, \code{\link{plot.leverage.ppm}},
  \code{\link{[.leverage.ppm}}, \code{\link{as.function.leverage.ppm}}.
}
\examples{
 fit <- ppm(cells ~ x)
 a <- leverage(fit)
 integral(a)
}
\author{
  \spatstatAuthors
}
\keyword{spatial}
\keyword{manip}
