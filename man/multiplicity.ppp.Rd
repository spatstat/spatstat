\name{multiplicity.ppp}
\alias{multiplicity}
\alias{multiplicity.default}
\alias{multiplicity.data.frame}
\alias{multiplicity.ppp}
\alias{multiplicity.ppx}
\title{Count Multiplicity of Duplicate Points}
\description{
  Counts the number of duplicates for each point
  in a spatial point pattern.
}
\usage{
 multiplicity(x)

 \method{multiplicity}{ppp}(x)

 \method{multiplicity}{ppx}(x)

 \method{multiplicity}{data.frame}(x)

 \method{multiplicity}{default}(x)
}
\arguments{
  \item{x}{
    A spatial point pattern
    (object of class \code{"ppp"} or \code{"ppx"})
    or a vector, matrix or data frame.
  }
}
\value{
  A vector of integers (multiplicities) of length equal to the
  number of points in \code{x}.
}
\details{
  Two points in a point pattern are deemed to be identical
  if their \eqn{x,y} coordinates are the same,
  and their marks are also the same (if they carry marks).
  The Examples section illustrates how it is possible for
  a point pattern to contain a pair of identical points.

  For each point in \code{x}, the function \code{multiplicity} counts how many
  points are identical to it, and returns the vector of counts.

  The argument \code{x} can also be a vector, a matrix or a data frame.
  When \code{x} is a vector, \code{m <- multiplicity(x)} is a vector
  of the same length as \code{x}, and \code{m[i]} is the
  number of elements of \code{x} that are identical to \code{x[i]}.
  When \code{x} is a matrix or data frame,
  \code{m <- multiplicity(x)} is a vector
  of length equal to the number of rows of \code{x}, and \code{m[i]}
  is the number of rows of \code{x} that are identical to
  the \code{i}th row.
}
\seealso{
  \code{\link{ppp.object}},
  \code{\link{duplicated.ppp}},
  \code{\link{unique.ppp}}
}
\examples{
   X <- ppp(c(1,1,0.5,1), c(2,2,1,2), window=square(3), check=FALSE)
   m <- multiplicity(X)

   # unique points in X, marked by their multiplicity
   first <- !duplicated(X)
   Y <- X[first] \%mark\% m[first]
}

\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/},
  Rolf Turner
  \email{r.turner@auckland.ac.nz}
  and Sebastian Meyer.
}
\keyword{spatial}
\keyword{utilities}
 
