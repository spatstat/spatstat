\name{methods.units}
\Rdversion{1.1}
\alias{methods.units} %DoNotExport
\alias{print.units}
\alias{summary.units}
\alias{rescale.units}
\alias{compatible.units}
\title{
  Methods for Units
}
\description{
  Methods for class \code{"units"}.
}
\usage{
  \method{print}{units}(x, ...)
  \method{summary}{units}(object, ...)
  \method{rescale}{units}(X, s, unitname)
  \method{compatible}{units}(A,B, ..., coerce=TRUE)
}
\arguments{
  \item{x,X,A,B,object}{
    Objects of class \code{"units"} representing
    units of length.
  }
  \item{s}{Conversion factor: the new units are \code{s} times the old units.}
  \item{\dots}{
    Other arguments.
    For \code{print.units} these arguments are passed to \code{print.default}.
    For \code{summary.units} they are ignored.
    For \code{compatible.units} these arguments
    are other objects of class \code{"units"}.
  }
  \item{coerce}{
    Logical. If \code{TRUE}, a null unit of length is compatible with
    any non-null unit.
  }
  \item{unitname}{Optional new name for the unit. If present, this overrides the
    rescaling operation and simply substitutes the new name for the old one.}
}
\details{
  These are methods for the generic functions
  \code{\link{print}}, \code{\link{summary}}, \code{\link{rescale}}
  and \code{\link{compatible}}
  for the class \code{"units"}.

  An object of class \code{"units"} represents a unit of length.

  The \code{print} method prints a description of the unit of length,
  and the \code{summary} method gives a more detailed description.

  The \code{rescale} method changes the unit of length by rescaling it.
  
  The \code{compatible} method tests whether two or more units of length
  are compatible.
}
\value{
  For \code{print.units} the value is \code{NULL}.
  For \code{summary.units} the value is an object of class
  \code{summary.units} (with its own print method).
  For \code{rescale.units} the value is another object of class \code{"units"}.
  For \code{compatible.units} the result is logical.
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\seealso{
  \code{\link{box3}},
  \code{\link{print}},
  \code{\link{unitname}}
}
\keyword{spatial}
\keyword{methods}
