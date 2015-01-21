\name{solapply}
\alias{solapply}
\alias{anylapply}
\title{
  Apply a Function Over a List and Obtain a List of Objects
}
\description{
  Applies the function \code{FUN} to each element of the list \code{X},
  and returns the result as a list of class \code{"solist"}
  or \code{"anylist"} as appropriate.
}
\usage{
  anylapply(X, FUN, \dots)

  solapply(X, FUN, \dots, check = TRUE, promote = TRUE, demote = FALSE)
}
\arguments{
  \item{X}{A list.}
  \item{FUN}{
    Function to be applied to each element of \code{X}.
  }
  \item{\dots}{
    Additional arguments to \code{FUN}.
  }
  \item{check,promote,demote}{
    Arguments passed to \code{\link{solist}} which determine
    how to handle different classes of objects.
  }
}
\details{
  These convenience functions are similar to \code{\link[base]{lapply}}
  except that they return a list of class \code{"solist"} or
  \code{"anylist"}.

  In both functions, the result is computed by 
  \code{lapply(X, FUN, \dots)}.

  In \code{anylapply} the result is converted to a list of class
  \code{"anylist"} and returned.

  In \code{solapply} the result is converted to
  a list of class \code{"solist"} \bold{if possible}, using 
  \code{\link{as.solist}}. If this is not possible,
  then the behaviour depends on the argument \code{demote}.
  If \code{demote=TRUE} the result will be returned as a
  list of class \code{"anylist"}. If \code{demote=FALSE} (the default),
  an error occurs.
}
\value{
  A list, usually of class \code{"solist"}.
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  Rolf Turner
  \email{r.turner@auckland.ac.nz}
  and Ege Rubak
  \email{rubak@math.aau.dk}
}
\seealso{
  \code{\link{solist}}, 
  \code{\link{anylist}}.
}
\examples{
  solapply(waterstriders, density)
}
\keyword{spatial}
\keyword{list}
\keyword{manip}
