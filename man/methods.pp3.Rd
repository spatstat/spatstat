\name{methods.pp3}
\Rdversion{1.1}
\alias{methods.pp3} %DoNotExport
\alias{print.pp3}
\alias{summary.pp3}
\alias{print.summary.pp3}
\alias{unitname.pp3}
\alias{unitname<-.pp3}
\title{
  Methods for three-dimensional point patterns
}
\description{
  Methods for class \code{"pp3"}.
}
\usage{
  \method{print}{pp3}(x, ...)
  \method{print}{summary.pp3}(x, ...)
  \method{summary}{pp3}(object, ...)
  \method{unitname}{pp3}(x)
  \method{unitname}{pp3}(x) <- value
}
\arguments{
  \item{x,object}{
    Object of class \code{"pp3"}.
  }
  \item{\dots}{
    Ignored.
  }
  \item{value}{
    Name of the unit of length. See \code{\link{unitname}}.
  }
}
\details{
  These are methods for the generic functions
  \code{\link{print}}, \code{\link{summary}},
  \code{\link{unitname}} and  \code{\link{unitname<-}}
  for the class \code{"pp3"} of three-dimensional point patterns.

  The \code{print} and \code{summary} methods print a description
  of the point pattern.

  The \code{unitname} method extracts the name of the unit of
  length in which the point coordinates are expressed.
  The \code{unitname<-} method assigns the name of the unit of length.
}
\value{
  For \code{print.pp3} the value is \code{NULL}.
  For \code{unitname.pp3} an object of class \code{"units"}.
}
\author{Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  and Rolf Turner \email{r.turner@auckland.ac.nz}
  
}
\seealso{
  \code{\link{pp3}},
  \code{\link{print}},
  \code{\link{unitname}}
  \code{\link{unitname<-}}
}
\examples{
   X <- pp3(runif(42),runif(42),runif(42), box3(c(0,1), unitname="mm"))
   X
   unitname(X)
   unitname(X) <- c("foot", "feet")
   summary(X)
}
\keyword{spatial}
\keyword{methods}

