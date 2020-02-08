\name{methods.boxx}
\Rdversion{1.1}
\alias{methods.boxx}  %DoNotExport
\alias{is.boxx}
\alias{print.boxx}
\alias{unitname.boxx}
\alias{unitname<-.boxx}
\title{
  Methods for Multi-Dimensional Box
}
\description{
  Methods for class \code{"boxx"}.
}
\usage{
  \method{is}{boxx}(x)
  \method{print}{boxx}(x, ...)
  \method{unitname}{boxx}(x)
  \method{unitname}{boxx}(x) <- value
}
\arguments{
  \item{x}{
    Object of class \code{"boxx"} representing a multi-dimensional box.
  }
  \item{\dots}{
    Other arguments passed to \code{print.default}.
  }
  \item{value}{
    Name of the unit of length. See \code{\link{unitname}}.
  }
}
\details{
  These are methods for the generic functions
  \code{\link{is}}, \code{\link{print}} and \code{\link{unitname}}
  for the class \code{"boxx"} of multi-dimensional boxes.

  The \code{is} method checks that the object is of class \code{"boxx"}, the
  \code{print} method prints a description of the box,
  the \code{unitname} method extracts the name of the unit of
  length in which the box coordinates are expressed, while the assignment method
  for \code{unitname} assigns this unit name.
}
\value{
  For \code{is.boxx} the value is a logical.
  For \code{print.boxx} the value is \code{NULL}.
  For \code{unitname.boxx} an object of class \code{"units"}.
  For \code{unitname<-.boxx} the result is the updated \code{"boxx"} object \code{x}.
}
\author{
  \spatstatAuthors
}
\seealso{
  \code{\link{boxx}},
  \code{\link{is}},
  \code{\link{print}},
  \code{\link{unitname}}
}
\examples{
   X <- boxx(c(0,10),c(0,10),c(0,5),c(0,1), unitname=c("metre", "metres"))
   X
   is.boxx(X)
   unitname(X)
   # Northern European usage
   unitname(X) <- "meter"
}
\keyword{spatial}
\keyword{methods}
