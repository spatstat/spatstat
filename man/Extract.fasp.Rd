\name{Extract.fasp}
\alias{[.fasp}
\title{Extract Subset of Function Array}
\description{
  Extract a subset of a function array (an object of class
  \code{"fasp"}).
}
\usage{
  \method{[}{fasp}(x, I, J, drop=TRUE,\dots)
}
\arguments{
  \item{x}{
    A function array. An object of class \code{"fasp"}.
  }
  \item{I}{
    any valid expression for a subset of the row indices of the array.
  }
  \item{J}{
    any valid expression for a subset of the column indices of the array.
  }
  \item{drop}{
      Logical.
      When the selected subset consists of only one cell of the array,
      if \code{drop=FALSE} the result is still returned as a
      \eqn{1 \times 1}{1 * 1} array of functions (class \code{"fasp"})
      while if \code{drop=TRUE} it is returned as a
      function (class \code{"fv"}).
  }
  \item{\dots}{Ignored.}
} 
\value{
  A function array (of class \code{"fasp"}).
  Exceptionally, if the array has only one cell, and
  if \code{drop=TRUE}, then the result is a function value table
  (class \code{"fv"}).
}
\details{
  A function array can be regarded as a matrix whose entries
  are functions. See \code{\link{fasp.object}} for an explanation of
  function arrays.

  This routine extracts a sub-array according to the usual
  conventions for matrix indexing.
}
\seealso{
  \code{\link{fasp.object}}
}
\examples{
 # Lansing woods data - multitype points with 6 types
 data(lansing)
 \testonly{
 # smaller dataset
   lansing <- lansing[ seq(1,lansing$n,by=45), ]
 }
 # compute 6 x 6 array of all cross-type K functions
 a <- alltypes(lansing, "K")

 # extract first three marks only
 b <- a[1:3,1:3]
 \dontrun{plot(b)}
 # subset pertaining to hickories
 h <- a[levels(lansing$marks) == "hickory", ]
 \dontrun{plot(h)}
}
\author{
  Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au},
 Rolf Turner \email{r.turner@auckland.ac.nz}
 and Ege Rubak \email{rubak@math.aau.dk}
}
\keyword{spatial}
\keyword{manip}
