\name{closetriples}
\alias{closetriples}
\title{
  Close Triples of Points
}
\description{
  Low-level function to find all close triples of points.
}
\usage{
closetriples(X, rmax)
}
\arguments{
  \item{X}{
    Point pattern (object of class \code{"ppp"} or \code{"pp3"}).
  }
  \item{rmax}{
    Maximum distance between each pair of points in a triple.
  }
}
\details{
  This low-level function
  finds all triples of points in a point pattern
  in which each pair lies closer than \code{rmax}.
}
\value{
  A data frame with columns \code{i,j,k} giving the indices of the
  points in each triple, and a column \code{diam} giving the diameter
  (maximum pairwise distance) in the triple.
}
\author{
  \spatstatAuthors.
}
\seealso{
  \code{\link{closepairs}},
  \code{\link{Tstat}}.
}
\examples{
   closetriples(redwoodfull, 0.02)
   closetriples(redwoodfull, 0.005)
}
\keyword{spatial}
\keyword{math}
