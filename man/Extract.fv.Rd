\name{Extract.fv}
\alias{[.fv}
\title{Extract Subset of Function Values}
\description{
  Extract a subset of an object of class \code{"fv"}.
}
\usage{
  \method{[}{fv}(x, i, j, \dots, drop=FALSE)
}
\arguments{
  \item{x}{
    a function value object,
    of class \code{"fv"} (see \code{\link{fv.object}}).
    Essentially a data frame.
  }
  \item{i}{
    any appropriate subset index.
    Selects a subset of the rows of the data frame, i.e.
    a subset of the domain of the function(s) represented by \code{x}.
  }
  \item{j}{
    any appropriate subset index for the columns of the data frame.
    Selects some of the functions present in \code{x}.
  }
  \item{drop, \dots}{Ignored.}
} 
\value{
  A function value object (of class \code{"fv"}).
}
\details{
  This is a method for \code{"["} for the class \code{"fv"}.
  It is very similar to \code{\link{[.data.frame}}
  except for a few extra checks on the sanity of the result.
}
\seealso{
  \code{\link{fv.object}}
}
\examples{
 data(cells)
 K <- Kest(cells)

 # discard the estimates of K(r) for r  > 0.1
 Ksub <- K[K$r <= 0.1, ]

 # discard the border method estimator
 Ksub <- K[ , names(K) != "border"]
 
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{manip}
