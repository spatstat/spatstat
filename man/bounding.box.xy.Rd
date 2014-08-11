\name{bounding.box.xy}
\alias{bounding.box.xy}
\title{Convex Hull of Points}
\description{
  Computes the smallest rectangle containing a set of points.
}
\usage{
bounding.box.xy(x, y=NULL)
}
\arguments{
  \item{x}{
    vector of \code{x} coordinates of observed points,
    or a 2-column matrix giving \code{x,y} coordinates,
    or a list with components \code{x,y} giving coordinates
    (such as a point pattern object of class \code{"ppp"}.)
  }
  \item{y}{(optional) vector of \code{y} coordinates of observed points,
    if \code{x} is a vector.}
}
\value{
  A window (an object of class \code{"owin"}).
}
\details{
  Given an observed pattern of points with coordinates 
  given by \code{x} and \code{y}, this function finds the smallest
  rectangle, with sides parallel to the coordinate axes, that contains
  all the points, and returns it as a window.
}
\seealso{
  \code{\link{owin}},
  \code{\link{as.owin}},
  \code{\link{convexhull.xy}},
  \code{\link{ripras}}
}
\examples{
  x <- runif(30)
  y <- runif(30)
  w <- bounding.box.xy(x,y)
  plot(owin(), main="bounding.box.xy(x,y)")
  plot(w, add=TRUE)
  points(x,y)

  X <- rpoispp(30)
  plot(X, main="bounding.box.xy(X)")
  plot(bounding.box.xy(X), add=TRUE)
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{utilities}




