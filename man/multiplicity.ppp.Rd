\name{multiplicity.ppp}
\alias{multiplicity}
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
}
\arguments{
  \item{x}{
    A spatial point pattern
    (object of class \code{"ppp"} or \code{"ppx"})
    or some other type of spatial object.
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

  For each point in \code{x}, this function counts how many
  points are identical to it, and returns the vector of counts.
}
\seealso{
  \code{\link{ppp.object}},
  \code{\link{duplicated.ppp}},
  \code{\link{unique.ppp}}
}
\examples{
   X <- ppp(c(1,1,0.5), c(2,2,1), window=square(3))
   m <- multiplicity(X)

   # unique points in X, marked by multiplicity
   first <- !duplicated(X)
   Y <- X[first] \%mark\% m[first]
}

\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{utilities}
 
