\name{ripras}
\alias{ripras}
\title{Estimate window from points alone}
\description{
  Given an observed pattern of points,
  computes the Ripley-Rasson estimate of 
  the spatial domain from which they came.
}
\usage{
 ripras(x, y=NULL, shape="convex", f)
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
  \item{shape}{String indicating the type of window to be estimated:
    either \code{"convex"} or \code{"rectangle"}.
  }
  \item{f}{
    (optional) scaling factor. See Details.
  }
}
\value{
  A window (an object of class \code{"owin"}).
}
\details{
  Given an observed pattern of points with coordinates 
  given by \code{x} and \code{y}, this function computes 
  an estimate due to Ripley and Rasson (1977) of the
  spatial domain from which the points came. 

  The points are
  assumed to have been generated independently and uniformly
  distributed inside an unknown domain \eqn{D}.

  If \code{shape="convex"} (the default), the domain \eqn{D} is assumed
  to be a convex set. The maximum
  likelihood estimate of \eqn{D} is the convex hull of the 
  points (computed by \code{\link{convexhull.xy}}).
  Analogously to the problems of estimating the endpoint
  of a uniform distribution, the MLE is not optimal.
  Ripley and Rasson's estimator is a rescaled copy of the convex hull,
  centred at the centroid of the convex hull.
  The scaling factor is 
  \eqn{1/sqrt(1 - m/n)}{1/\sqrt{1 - \frac m n}}
  where \eqn{n} is the number of data points and 
  \eqn{m} the number of vertices of the convex hull.
  The scaling factor may be overridden using the argument \code{f}.

  If \code{shape="rectangle"}, the domain \eqn{D} is assumed
  to be a rectangle with sides parallel to the coordinate axes. The maximum
  likelihood estimate of \eqn{D} is the bounding box of the points
  (computed by \code{\link{bounding.box.xy}}). The Ripley-Rasson
  estimator is a rescaled copy of the bounding box,
  with scaling factor \eqn{(n+1)/(n-1)}
  where \eqn{n} is the number of data points,
  centred at the centroid of the bounding box.
  The scaling factor may be overridden using the argument \code{f}.
}
\seealso{
  \code{\link{owin}},
  \code{\link{as.owin}},
  \code{\link{bounding.box.xy}},
  \code{\link{convexhull.xy}}
}
\examples{
  x <- runif(30)
  y <- runif(30)
  w <- ripras(x,y)
  plot(owin(), main="ripras(x,y)")
  plot(w, add=TRUE)
  points(x,y)

  X <- rpoispp(15)
  plot(X, main="ripras(X)")
  plot(ripras(X), add=TRUE)

  # two points insufficient
  ripras(c(0,1),c(0,0))
  # triangle
  ripras(c(0,1,0.5), c(0,0,1))
  # three collinear points
  ripras(c(0,0,0), c(0,1,2))
}
\references{
  Ripley, B.D. and Rasson, J.-P. (1977)
  Finding the edge of a Poisson forest.
  \emph{Journal of Applied Probability},
  \bold{14}, 483 -- 491.
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{utilities}




