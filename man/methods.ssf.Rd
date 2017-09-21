\name{methods.ssf}
\alias{methods.ssf} %DoNotExport
\alias{marks.ssf}
\alias{marks<-.ssf}
\alias{unmark.ssf}
\alias{as.im.ssf}
\alias{as.function.ssf}
\alias{as.ppp.ssf}
\alias{print.ssf}
\alias{range.ssf}
\alias{min.ssf}
\alias{max.ssf}
\alias{integral.ssf}
\title{Methods for Spatially Sampled Functions}
\description{
  Methods for various generic commands, for the class
  \code{"ssf"} of spatially sampled functions.
}
\usage{
  \method{marks}{ssf}(x, \dots)

  \method{marks}{ssf}(x, \dots) <- value

  \method{unmark}{ssf}(X)

  \method{as.im}{ssf}(X, \dots)

  \method{as.function}{ssf}(x, \dots)

  \method{as.ppp}{ssf}(X, \dots)

  \method{print}{ssf}(x, \dots, brief=FALSE)

  \method{range}{ssf}(x, \dots)

  \method{min}{ssf}(x, \dots)

  \method{max}{ssf}(x, \dots)

  \method{integral}{ssf}(f, domain=NULL, ..., weights=attr(f, "weights"))
}
\arguments{
  \item{x,X,f}{
    A spatially sampled function (object of class \code{"ssf"}). 
  }
  \item{\dots}{Arguments passed to the default method.}
  \item{brief}{Logical value controlling the amount of detail printed.}
  \item{value}{Matrix of replacement values for the function.}
  \item{domain}{Optional.
    Domain of integration. An object of class\code{"owin"}.
  }
  \item{weights}{
    Optional. Numeric vector of weights associated with the
    sample points.
  }
}
\value{
  \code{marks} returns a matrix.

  \code{marks(x) <- value} returns an object of class \code{"ssf"}.
  
  \code{as.owin} returns a window (object of class \code{"owin"}).

  \code{as.ppp} and \code{unmark}
  return a point pattern (object of class \code{"ppp"}).

  \code{as.function} returns a \code{function(x,y)} of class \code{"funxy"}.
  
  \code{print} returns \code{NULL}.

  \code{range} returns a numeric vector of length 2.
  \code{min} and \code{max} return a single numeric value.

  \code{integral} returns a numeric value (if \code{x} had numeric
  values) or a numeric vector (if \code{x} had vector values).
}
\details{
  An object of class \code{"ssf"} represents a
  function (real- or vector-valued) that has been
  sampled at a finite set of points.

  The commands documented here are methods for this class,
  for the generic commands
  \code{\link[spatstat]{marks}},
  \code{\link[spatstat]{marks<-}},
  \code{\link[spatstat]{unmark}},
  \code{\link[spatstat]{as.im}},
  \code{\link{as.function}},
  \code{\link[spatstat]{as.ppp}},
  \code{\link{print}},
  \code{\link{range}},
  \code{\link{min}},
  \code{\link{max}}
  and \code{\link[spatstat]{integral}}.
}
\seealso{
  \code{\link{ssf}}
}
\examples{
  X <- cells[1:4]
  f <- ssf(X, nndist(X, k=1:3))
  f
  marks(f)
  as.ppp(f)
  as.im(f)
}
\author{Adrian Baddeley}
\keyword{spatial}
\keyword{methods}

