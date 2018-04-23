\name{methods.distfun}
\alias{methods.distfun} %DoNotExport
\alias{shift.distfun}
\alias{rotate.distfun}
\alias{scalardilate.distfun}
\alias{affine.distfun}
\alias{flipxy.distfun}
\alias{reflect.distfun}
\alias{rescale.distfun}
\Rdversion{1.1}
\title{
  Geometrical Operations for Distance Functions
}
\description{
  Methods for objects of the class \code{"distfun"}.
}
\usage{
\method{shift}{distfun}(X, \dots)

\method{rotate}{distfun}(X, \dots)

\method{scalardilate}{distfun}(X, \dots)

\method{affine}{distfun}(X, \dots)

\method{flipxy}{distfun}(X)

\method{reflect}{distfun}(X)

\method{rescale}{distfun}(X, s, unitname)
}
\arguments{
  \item{X}{
    Object of class \code{"distfun"} representing the distance function
    of a spatial object.
  }
  \item{\dots}{
    Arguments passed to the next method for the geometrical operation.
    See Details.
  }
  \item{s, unitname}{
    Arguments passed to the next method for \code{\link{rescale}}.
  }
}
\details{
  These are methods for the generic functions
  \code{\link{shift}},
  \code{\link{rotate}},
  \code{\link{scalardilate}},
  \code{\link{affine}},
  \code{\link{flipxy}} and \code{\link{reflect}}
  which perform geometrical operations on spatial objects,
  and for the generic \code{\link{rescale}} which changes the unit of length.

  The argument \code{X} should be an object of class \code{"distfun"}
  representing the distance function of a spatial object \code{Y}.
  Objects of class \code{"distfun"} are created by \code{\link{distfun}}.

  The methods apply the specified geometrical transformation to
  the original object \code{Y}, producing a new object \code{Z}
  of the same type as \code{Y}. They then create a new \code{distfun}
  object representing the distance function of \code{Z}.
}
\value{
  Another object of class \code{"distfun"}.
}
\author{
  \spatstatAuthors.
}
\seealso{
  \code{\link{distfun}},
  \code{\link{methods.funxy}}.
}
\examples{
   (f <- distfun(letterR))
   plot(f)
   flipxy(f)
   shift(f, origin="midpoint")
   plot(rotate(f, angle=pi/2))

   (g <- distfun(lansing))
   rescale(g)
}
\keyword{spatial}
\keyword{methods}
