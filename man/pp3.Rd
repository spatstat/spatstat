\name{pp3}
\Rdversion{1.1}
\alias{pp3}
\title{
  Three Dimensional Point Pattern
}
\description{
  Create a three-dimensional point pattern
}
\usage{
pp3(x, y, z, ...)
}
\arguments{
  \item{x,y,z}{
    Numeric vectors of equal length, containing Cartesian coordinates
    of points in three-dimensional space.
}
  \item{\dots}{
    Arguments passed to \code{\link{as.box3}} to determine
    the three-dimensional box in which the points have been observed.
  }
}
\details{
  An object of class \code{"pp3"} represents a pattern of points in
  three-dimensional space. The points are assumed to have been observed
  by exhaustively inspecting a three-dimensional rectangular box. The
  boundaries of the box are included as part of the dataset.
}
\value{
  Object of class \code{"pp3"} representing a three dimensional
  point pattern. Also belongs to class \code{"ppx"}. 
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\seealso{
  \code{\link{box3}}, 
  \code{\link{print.pp3}}, 
  \code{\link{ppx}} 
}
\examples{
   X <- pp3(runif(10), runif(10), runif(10), box3(c(0,1)))
}
\keyword{spatial}
\keyword{datagen}
