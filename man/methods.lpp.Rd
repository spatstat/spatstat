\name{methods.lpp}
\alias{methods.lpp} %DoNotExport
\Rdversion{1.1}
\alias{as.ppp.lpp}
\alias{as.psp.lpp}
\alias{marks<-.lpp}
\alias{nsegments.lpp}
\alias{print.lpp}
\alias{print.summary.lpp}
\alias{summary.lpp}
\alias{unitname.lpp}
\alias{unitname<-.lpp}
\alias{unmark.lpp}
\title{
  Methods for Point Patterns on a Linear Network
}
\description{
  These are methods specifically for the class \code{"lpp"} of point patterns on
  linear networks.
}
\usage{
\method{as.ppp}{lpp}(X, ..., fatal=TRUE)

\method{as.psp}{lpp}(x, ..., fatal=TRUE)

\method{marks}{lpp}(x, ...) <- value

\method{nsegments}{lpp}(x)

\method{print}{lpp}(x, ...)

\method{print}{summary.lpp}(x, ...)

\method{summary}{lpp}(object, ...)

\method{unitname}{lpp}(x)

\method{unitname}{lpp}(x) <- value

\method{unmark}{lpp}(X)
}
\arguments{
  \item{x,X,object}{
    An object of class \code{"lpp"} representing a point pattern
    on a linear network.
  }
  \item{\dots}{
    Arguments passed to other methods.
  }
  \item{value}{
    Replacement value for the \code{marks} or \code{unitname}
    of \code{x}. See Details.
  }
  \item{fatal}{
    Logical value indicating whether data in the wrong format
    should lead to an error (\code{fatal=TRUE}) or a warning
    (\code{fatal=FALSE}).
  }
}
\details{
  These are methods for the generic functions
  \code{\link{as.ppp}},
  \code{\link{as.psp}},
  \code{\link{marks<-}},
  \code{\link{nsegments}},
  \code{\link{print}},
  \code{\link{summary}},
  \code{\link{unitname}}, 
  \code{\link{unitname<-}} and
  \code{\link{unmark}}
  for objects of the class \code{"lpp"}.

  For \code{"marks<-.lpp"} the replacement \code{value}
  should be either \code{NULL}, or a vector of length equal
  to the number of points in \code{x},
  or a data frame with one row for each point in \code{x}.
  
  For \code{"unitname<-.lpp"} the replacement \code{value}
  should be a valid name for the unit of length, as
  described in \code{\link{unitname}}.
}
\section{Other methods}{
  An object of class \code{"lpp"} also inherits the class
  \code{"ppx"} for which many other methods are available.
  See \code{\link[spatstat:methods.ppx]{methods.ppx}}.
}
\value{
  See the documentation on the corresponding generic function.
}
\author{
  \adrian
}
\seealso{
  \code{\link{lpp}},
  \code{\link{intensity.lpp}},
  \code{\link[spatstat:methods.ppx]{methods.ppx}}
}
\examples{
  X <- runiflpp(10, simplenet)
  X
  as.ppp(X)
  summary(X)
  unitname(X) <- c("furlong", "furlongs")
}
\keyword{spatial}
\keyword{methods}
