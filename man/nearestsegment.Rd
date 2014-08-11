\name{nearestsegment}
\alias{nearestsegment}
\title{Find Line Segment Nearest to Each Point}
\description{
  Given a point pattern and a line segment pattern,
  this function finds the nearest line segment for each point.
}
\usage{
nearestsegment(X, Y)
}
\arguments{
  \item{X}{A point pattern (object of class \code{"ppp"}).}
  \item{Y}{A line segment pattern (object of class \code{"psp"}).}
}
\details{
  The distance between a point \code{x} and a straight line segment \code{y} is
  defined to be the shortest Euclidean distance between \code{x} and any
  location on \code{y}. This algorithm first calculates the distance
  from each point of \code{X} to
  each segment of \code{Y}. Then it determines, for each point \code{x} in
  \code{X}, which segment of \code{Y} is closest. The index of this
  segment is returned.
}
\value{
  Integer vector \code{v} (of length equal to the number of points in
  \code{X}) identifying the nearest segment to each point.
  If \code{v[i] = j}, then
  \code{Y[j]} is the line segment lying closest to \code{X[i]}.
}
\author{
  Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\seealso{
  \code{\link{project2segment}} to project each point of \code{X} to
  a point lying on one of the line segments.

  Use \code{\link{distmap.psp}} to 
  identify the nearest line segment for each pixel in a grid.
}
\examples{
  X <- runifpoint(3)
  Y <- as.psp(matrix(runif(20), 5, 4), window=owin())
  v <- nearestsegment(X,Y)
  plot(Y)
  plot(X, add=TRUE)
  plot(X[1], add=TRUE, col="red")
  plot(Y[v[1]], add=TRUE, lwd=2, col="red") 
}
\keyword{spatial}
\keyword{math}
