\name{inside.owin}
\alias{inside.owin}
\title{Test Whether Points Are Inside A Window}
\description{
 Test whether points lie inside or outside
 a given window.
}
\usage{
 inside.owin(x, y, w)
}
\arguments{
  \item{x}{
    Vector of \eqn{x} coordinates of points to be tested.
    (Alternatively, a point pattern object providing both
    \eqn{x} and \eqn{y} coordinates.)
  }
  \item{y}{
    Vector of \eqn{y} coordinates of points to be tested.
  }
  \item{w}{A window.
    This should be an object of class \code{\link{owin}},
    or can be given in any format acceptable to \code{\link{as.owin}()}.
  }
}
\value{
  Logical vector whose \code{i}th entry is 
  \code{TRUE} if the corresponding point \code{(x[i],y[i])}
  is inside \code{w}. 
}
\details{
  This function tests whether each of the points 
  \code{(x[i],y[i])} lies inside or outside
  the window \code{w} and returns \code{TRUE} if it is inside.
  
  The boundary of the window is treated as being inside.

  If \code{w} is of type \code{"rectangle"} or 
  \code{"polygonal"}, the algorithm uses analytic geometry
  (the discrete Stokes theorem).
  Computation time is linear in the number of points
  and (for polygonal windows) in the number of vertices of the
  boundary polygon. Boundary cases are correct to single
  precision accuracy.
  
  If \code{w} is of type \code{"mask"} then the 
  pixel closest to \code{(x[i],y[i])} is tested. The
  results may be incorrect for points lying within 
  one pixel diameter of the window boundary.

  Normally \code{x} and \code{y} must be numeric vectors of
  equal length (length zero is allowed) containing the coordinates
  of points. Alternatively \code{x}
  can be a point pattern (object of class \code{"ppp"}) while \code{y}
  is missing; then the coordinates of the point pattern are extracted.
}
\seealso{
  \code{\link{owin.object}},
  \code{\link{as.owin}}
}
\examples{
  # hexagonal window
  k <- 6
  theta <- 2 * pi * (0:(k-1))/k
  co <- cos(theta)
  si <- sin(theta)
  mas <- owin(c(-1,1), c(-1,1), poly=list(x=co, y=si))
  \dontrun{
  plot(mas)
  }

  # random points in rectangle
  x <- runif(30,min=-1, max=1)
  y <- runif(30,min=-1, max=1)

  ok <- inside.owin(x, y, mas)

  \dontrun{
  points(x[ok], y[ok])
  points(x[!ok], y[!ok], pch="x")
  }  
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{math}
