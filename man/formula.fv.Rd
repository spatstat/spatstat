\name{formula.fv}
\alias{formula.fv}
\alias{formula<-}
\alias{formula<-.fv}
\title{
  Extract or Change the Plot Formula for a Function Value Table
}
\description{
  Extract or change the default plotting formula
  for an object of class \code{"fv"} (function value table).
}
\usage{
\method{formula}{fv}(x, \dots)

formula(x, \dots) <- value

\method{formula}{fv}(x, \dots) <- value
}
\arguments{
  \item{x}{
    An object of class \code{"fv"},
    containing the values of several estimates of a function.
  }
  \item{\dots}{
    Arguments passed to other methods.
  }
  \item{value}{
    New value of the formula. Either a \code{formula} or a character
    string.
  }
}
\details{
  A function value table
  (object of class \code{"fv"}, see \code{\link{fv.object}})
  is a convenient way of storing and plotting
  several different estimates of the same function.

  The default behaviour of \code{plot(x)} for a function value table
  \code{x} is determined by a formula
  associated with \code{x} called its \emph{plot formula}.
  See \code{\link{plot.fv}} for explanation about these formulae.

  The function \code{formula.fv} is a method for the generic command
  \code{\link{formula}}. It extracts the plot formula associated with
  the object.

  The function \code{formula<-} is generic. It changes the formula
  associated with an object. 

  The function \code{formula<-.fv} is the method for \code{formula<-}
  for the class \code{"fv"}. It changes the plot formula associated with
  the object.
}
\value{
  The result of \code{formula.fv} is a character string containing the
  plot formula. The result of \code{formula<-.fv} is a new object of
  class \code{"fv"}.
}
\author{
  Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  and Rolf Turner \email{r.turner@auckland.ac.nz}
  
}
\seealso{
  \code{\link{fv}},  
  \code{\link{plot.fv}},
  \code{\link[stats]{formula}}.
}
\examples{
   K <- Kest(cells)
   formula(K)
   formula(K) <- (iso ~ r)
}
\keyword{spatial}
\keyword{methods}
