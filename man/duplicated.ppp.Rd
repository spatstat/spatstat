\name{duplicated.ppp}
\alias{duplicated.ppp}
\title{Determine Duplicated Points in a Spatial Point Pattern}
\description{
  Determines which points in a spatial point pattern
  are duplicates of previous points, and returns a logical vector.
}
\usage{
 \method{duplicated}{ppp}(x, \dots, rule=c("spatstat", "deldir"))
}
\arguments{
  \item{x}{
    A spatial point pattern
    (object of class \code{"ppp"}).
  }
  \item{\dots}{
    Ignored.
  }
  \item{rule}{
    Character string.
    The rule for determining duplicated points.
  }
}
\value{
  A logical vector of length equal to the number of points in \code{x}.
}
\details{
  This is a method for the generic function \code{duplicated} for 
  point pattern datasets (of class \code{"ppp"}, see \code{\link{ppp.object}}).

  Two points in a point pattern are deemed to be identical
  if their \eqn{x,y} coordinates are the same,
  and their marks are also the same (if they carry marks).
  The Examples section illustrates how it is possible for
  a point pattern to contain a pair of identical points.

  This function determines which points in \code{x} duplicate
  other points that appeared earlier in the sequence. It 
  returns a logical vector with entries that are \code{TRUE}
  for duplicated points and \code{FALSE} for unique (non-duplicated)
  points.

  If \code{rule="spatstat"} (the default), duplicated points are
  determined by testing equality of their coordinates and marks
  using \code{==}. This is the most stringent possible test.
  If \code{rule="deldir"}, duplicated points are
  determined using the function \code{\link[deldir]{duplicatedxy}}
  in the package \pkg{deldir}, which currently uses
  \code{\link{duplicated.data.frame}}. Setting \code{rule="deldir"}
  will ensure consistency with functions in the \pkg{deldir} package.
}
\seealso{
  \code{\link{ppp.object}},
  \code{\link{unique.ppp}},
  \code{\link{multiplicity.ppp}}
}
\examples{
   X <- ppp(c(1,1,0.5), c(2,2,1), window=square(3))
   duplicated(X)
   duplicated(X, rule="deldir")
}

\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{methods}
 
