\name{convexhull.xy}
\alias{convexhull.xy}
\title{Convex Hull of Points}
\description{
  Computes the convex hull of a set of points in two dimensions.
}
\usage{
convexhull.xy(x, y=NULL)
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
  given by \code{x} and \code{y}, this function computes 
  the convex hull of the points, and returns it as a window.
}
\seealso{
  \code{\link{owin}},
  \code{\link{as.owin}},
  \code{\link{convexhull}},
  \code{\link{bounding.box.xy}},
  \code{\link{ripras}}
}
\examples{
  x <- runif(30)
  y <- runif(30)
  w <- convexhull.xy(x,y)
  plot(owin(), main="convexhull.xy(x,y)", lty=2)
  plot(w, add=TRUE)
  points(x,y)

  X <- rpoispp(30)
  plot(X, main="convexhull.xy(X)")
  plot(convexhull.xy(X), add=TRUE)
}
\author{Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  and Rolf Turner \email{r.turner@auckland.ac.nz}
  
}
\keyword{spatial}
\keyword{utilities}




