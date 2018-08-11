\name{methods.unitname}
\Rdversion{1.1}
\alias{methods.unitname} %DoNotExport
\alias{print.unitname}
\alias{summary.unitname}
\alias{rescale.unitname}
\alias{compatible.unitname}
\alias{harmonise.unitname}
\alias{harmonize.unitname}
\title{
  Methods for Units
}
\description{
  Methods for class \code{"unitname"}.
}
\usage{
  \method{print}{unitname}(x, ...)
  \method{summary}{unitname}(object, ...)
  \method{rescale}{unitname}(X, s, unitname)
  \method{compatible}{unitname}(A,B, ..., coerce=TRUE)
  \method{harmonise}{unitname}(..., coerce=TRUE, single=FALSE)
  \method{harmonize}{unitname}(..., coerce=TRUE, single=FALSE)
}
\arguments{
  \item{x,X,A,B,object}{
    Objects of class \code{"unitname"} representing
    units of length.
  }
  \item{\dots}{
    Other arguments.
    For \code{print.unitname} these arguments are passed to
    \code{\link[base]{print.default}}.
    For \code{summary.unitname} they are ignored.
    For \code{compatible.unitname} and \code{harmonise.unitname} these arguments
    are other objects of class \code{"unitname"}.
  }
  \item{s}{
    Conversion factor: the new units are \code{s} times the old units.
  }
  \item{unitname}{
    Optional new name for the unit. If present, this overrides the
    rescaling operation and simply substitutes the new name for the old
    one.
  }
  \item{coerce}{
    Logical. If \code{TRUE}, a null unit of length is compatible with
    any non-null unit.
  }
  \item{single}{
    Logical value indicating whether to return a single unitname,
    or a list of unitnames.
  }
}
\details{
  These are methods for the generic functions
  \code{\link{print}}, \code{\link{summary}}, \code{\link{rescale}}
  and \code{\link{compatible}}
  for the class \code{"unitname"}.

  An object of class \code{"unitname"} represents a unit of length.

  The \code{print} method prints a description of the unit of length,
  and the \code{summary} method gives a more detailed description.

  The \code{rescale} method changes the unit of length by rescaling it.
  
  The \code{compatible} method tests whether two or more units of length
  are compatible.

  The \code{harmonise} method returns the common unit of length
  if there is one. For consistency with other methods for
  \code{\link{harmonise}}, the result is a list of unitname objects,
  with one entry for each argument in \code{\dots}. All of these
  entries are identical. This can be overridden by setting
  \code{single=TRUE} when the result will be a single unitname object.
}
\value{
  For \code{print.unitname} the value is \code{NULL}.
  For \code{summary.unitname} the value is an object of class
  \code{summary.unitname} (with its own print method).
  For \code{rescale.unitname} the value is another object of class \code{"unitname"}.
  For \code{compatible.unitname} the result is logical.
  For \code{harmonise.unitname} the result is a list of
  identical unitnames if \code{single=FALSE} (the default),
  or a single unitname if \code{single=TRUE}.
}
\author{
  \spatstatAuthors.
}
\seealso{
  \code{\link{box3}},
  \code{\link{print}},
  \code{\link{unitname}}
}
\keyword{spatial}
\keyword{methods}
