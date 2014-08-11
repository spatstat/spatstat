\name{methods.linnet}
\alias{methods.linnet} %DoNotExport
\Rdversion{1.1}
\alias{as.linnet}
\alias{as.linnet.linnet}
\alias{as.linnet.lpp}
\alias{as.owin.linnet}
\alias{as.psp.linnet}
\alias{print.linnet}
\alias{summary.linnet}
\alias{unitname.linnet}
\alias{unitname<-.linnet}
\title{
  Methods for Linear Networks
}
\description{
  These are methods for the class \code{"linnet"} of linear networks.
}
\usage{
as.linnet(X, ...)
\method{as.linnet}{linnet}(X, ...)
\method{as.linnet}{lpp}(X, ..., fatal=TRUE)
\method{as.owin}{linnet}(W, ...)
\method{as.psp}{linnet}(x, ..., fatal=TRUE)
\method{print}{linnet}(x, ...)
\method{summary}{linnet}(object, ...)
\method{unitname}{linnet}(x)
\method{unitname}{linnet}(x) <- value
}
\arguments{
  \item{x,X,object,W}{
    An object of class \code{"linnet"} representing 
    a linear network.
  }
  \item{\dots}{
    Arguments passed to other methods.
  }
  \item{value}{
    A valid name for the unit of length for \code{x}.
    See \code{\link{unitname}}.
  }
  \item{fatal}{
    Logical value indicating whether data in the wrong format
    should lead to an error (\code{fatal=TRUE}) or a warning
    (\code{fatal=FALSE}).
  }
}
\details{
  The function \code{as.linnet} is generic.
  It converts data from some other format
  into an object of class \code{"linnet"}.
  The method \code{as.linnet.lpp} extracts the linear network
  information from an \code{lpp} object.
  
  The other functions are methods for the generic commands
  \code{\link{as.owin}},
  \code{\link{as.psp}},
  \code{\link{print}},
  \code{\link{summary}},
  \code{\link{unitname}} and
  \code{\link{unitname<-}}
  for the class \code{"linnet"}.

  The method \code{as.owin.linnet} extracts the window containing
  the linear network, and returns it as an object of class \code{"owin"}.
  
  The method \code{as.psp.linnet} extracts the
  lines of the linear network as a line segment pattern (object of class
  \code{"psp"}).
}
\value{
  For \code{as.linnet} the value is an object of class \code{"linnet"}.
  For other functions, see the help file for the corresponding
  generic function.
}
\author{
  Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
}
\seealso{
  \code{\link{linnet}}
}
\examples{
  data(simplenet)
  simplenet
  summary(simplenet)
  unitname(simplenet) <- c("cubit", "cubits")
}
\keyword{spatial}
\keyword{methods}
