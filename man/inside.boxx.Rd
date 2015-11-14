\name{inside.boxx}
\alias{inside.boxx}
\title{Test Whether Points Are Inside A Multidimensional Box}
\description{
 Test whether points lie inside or outside
 a given multidimensional box.
}
\usage{
 inside.boxx(\dots, w)
}
\arguments{
  \item{\dots}{
    Coordinates of points to be tested.
    One vector for each dimension (all of same length).
    (Alternatively, a single point pattern object
    of class \code{"\link{ppx}"}
    or its coordinates as a \code{"\link{hyperframe}"})
  }
  \item{w}{A window.
    This should be an object of class \code{\link{boxx}},
    or can be given in any format acceptable to \code{\link{as.boxx}()}.
  }
}
\value{
  Logical vector whose \code{i}th entry is 
  \code{TRUE} if the corresponding point is inside \code{w}. 
}
\details{
  This function tests whether each of the points 
  \code{(x[i],y[i])} lies inside or outside
  the window \code{w} and returns \code{TRUE} if it is inside.
  
  The boundary of the window is treated as being inside.

  Normally each argument provided (except \code{w}) must be numeric vectors of
  equal length (length zero is allowed) containing the coordinates
  of points.
  Alternatively a single point pattern (object of class \code{"ppx"})
  can be given; then the coordinates of the point pattern are extracted.
}
\seealso{
  \code{\link{boxx}},
  \code{\link{as.boxx}}
}
\examples{
  # Random points in box with side [0,2]
  w <- boxx(c(0,2), c(0,2), c(0,2))

  # Random points in box with side [-1,3]
  x <- runif(30, min=-1, max=3)
  y <- runif(30, min=-1, max=3)
  z <- runif(30, min=-1, max=3)

  # Points falling in smaller box
  ok <- inside.boxx(x, y, z, w=w)

  # Same using a point pattern as argument:
  X <- ppx(data = cbind(x, y, z), domain = boxx(c(0,3), c(0,3), c(0,3)))
  ok2 <- inside.boxx(X, w=w)
}
\author{
  \adrian
  
  
  \rolf
  
  and \ege
  
}
\keyword{spatial}
\keyword{math}
