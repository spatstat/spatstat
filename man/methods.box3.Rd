\name{methods.box3}
\Rdversion{1.1}
\alias{methods.box3}  %DoNotExport
\alias{print.box3}
\alias{unitname.box3}
\alias{unitname<-.box3}
\title{
  Methods for Three-Dimensional Box
}
\description{
  Methods for class \code{"box3"}.
}
\usage{
  \method{print}{box3}(x, ...)
  \method{unitname}{box3}(x) 
  \method{unitname}{box3}(x) <- value
}
\arguments{
  \item{x}{
    Object of class \code{"box3"} representing a three-dimensional box.
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
  \code{\link{print}} and \code{\link{unitname}}
  for the class \code{"box3"} of three-dimensional boxes.

  The \code{print} method prints a description of the box,
  while the \code{unitname} method extracts the name of the unit of
  length in which the box coordinates are expressed.
}
\value{
  For \code{print.box3} the value is \code{NULL}.
  For \code{unitname.box3} an object of class \code{"units"}.
}
\author{\adrian
  
  
  and \rolf
  
}
\seealso{
  \code{\link{box3}},
  \code{\link{print}},
  \code{\link{unitname}}
}
\examples{
   X <- box3(c(0,10),c(0,10),c(0,5), unitname=c("metre", "metres"))
   X
   unitname(X)
   # Northern European usage
   unitname(X) <- "meter"
}
\keyword{spatial}
\keyword{methods}
