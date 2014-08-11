\name{as.function.fv}
\alias{as.function.fv}
\alias{as.function.rhohat}
\title{
  Convert Function Value Table to Function
}
\description{
  Converts an object of class \code{"fv"} to an \R language function.
}
\usage{
  \method{as.function}{fv}(x, ..., value, extrapolate=FALSE)

  \method{as.function}{rhohat}(x, ..., value, extrapolate=TRUE)
}
\arguments{
  \item{x}{
    Object of class \code{"fv"} or \code{"rhohat"}.
}
  \item{\dots}{
    Ignored.
}
  \item{value}{
    Optional. String selecting one of the columns of \code{x}
    for use as the function value.
  }
  \item{extrapolate}{
    Logical value indicating whether to extrapolate the function
    outside the domain of \code{x}. See Details.
  }
}
\details{
  A function value table (object of class \code{"fv"})
  is a convenient way of storing and plotting
  several different estimates of the same function.
  Objects of this class are returned 
  by many commands in \pkg{spatstat}, such as \code{\link{Kest}}
  which returns an estimate of Ripley's \eqn{K}-function
  for a point pattern dataset.

  Sometimes it is useful to convert the function value table
  to a \code{function} in the \R language. This is done by
  \code{as.function.fv}. It converts an object \code{x} of class \code{"fv"}
  to an \R function \code{f}.

  If \code{f <- as.function(x)} then \code{f} is a function with one numeric
  argument, that performs linear interpolation between the values
  in the table \code{x}.

  Argument values lying outside the range of the table
  yield an \code{NA} value (if \code{extrapolate=FALSE})
  or the function value at the nearest endpoint of the range
  (if \code{extrapolate = TRUE}).

  The command \code{as.function.fv} is a method for the generic command
  \code{\link{as.function}}. 
}
\value{
  A \code{function} with one argument.
}
\author{
  Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\seealso{
   \code{\link{fv}}, 
   \code{\link{fv.object}}, 
   \code{\link{plot.fv}}, 
   \code{\link{Kest}}
}
\examples{
  data(cells)
  K <- Kest(cells)
  f <- as.function(K)
  f(0.1)
}
\keyword{spatial}
\keyword{methods}

