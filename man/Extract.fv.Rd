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
  \item{\dots}{
    Ignored.
  }
  \item{drop}{
     Logical. If \code{TRUE}, the result is a data frame or vector
     containing the selected rows and columns of data.
     If \code{FALSE} (the default),
     the result is another object of class \code{"fv"}.
  }
} 
\value{
  If \code{drop=FALSE}, a function value object (of class \code{"fv"}).
  If \code{drop=TRUE}, a data frame or vector.
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
 K <- Kest(cells)

 # discard the estimates of K(r) for r  > 0.1
 Ksub <- K[K$r <= 0.1, ]

 # discard the border method estimator
 Ksub <- K[ , names(K) != "border"]

 # read some values
 K[5, ,drop=TRUE]
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{manip}
