\name{stratrand}
\alias{stratrand}
\title{Stratified random point pattern}
\description{
  Generates a \dQuote{stratified random} pattern of points in a window,
  by dividing the window into rectangular tiles and placing
  \code{k} random points in each tile.
}
\usage{
 stratrand(window, nx, ny, k = 1)
}
\arguments{
  \item{window}{A window. 
    An object of class \code{\link{owin}},
    or data in any format acceptable to \code{\link{as.owin}()}.
  }
  \item{nx}{Number of tiles in each row.
  }
  \item{ny}{Number of tiles in each column.
  }
  \item{k}{Number of random points to generate in each tile.
  }
}
\value{
  A list with two components \code{x} and \code{y}, which are numeric
  vectors giving the coordinates of the random points.
}
\details{
  The bounding rectangle of \code{window} is divided into
  a regular \eqn{nx \times ny}{nx * ny} grid of rectangular tiles.
  In each tile, \code{k} random points are generated independently
  with a uniform distribution in that tile. 

  Note that some of these grid points may lie outside the window,
  if \code{window} is not of type \code{"rectangle"}. The function
  \code{\link{inside.owin}} can be used to select those grid points
  which do lie inside the window. See the examples.

  This function is useful in creating dummy points for quadrature
  schemes (see \code{\link{quadscheme}}) as well as in simulating
  random point patterns.
}
\seealso{
  \code{\link{quad.object}},
  \code{\link{quadscheme}},
  \code{\link{inside.owin}},
  \code{\link{gridcentres}}
}
\examples{
  w <- unit.square()
  xy <- stratrand(w, 10, 10)
  \dontrun{
  plot(w)
  points(xy)
  }

  # polygonal boundary
  bdry <- list(x=c(0.1,0.3,0.7,0.4,0.2),
               y=c(0.1,0.1,0.5,0.7,0.3))
  w <- owin(c(0,1), c(0,1), poly=bdry)
  xy <- stratrand(w, 10, 10, 3)
  \dontrun{
  plot(w)
  points(xy)
  }
  # determine which grid points are inside polygon
  ok <- inside.owin(xy$x, xy$y, w)
  \dontrun{
  plot(w)
  points(xy$x[ok], xy$y[ok])
  }
}
\author{Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  and Rolf Turner \email{r.turner@auckland.ac.nz}
  
}
\keyword{spatial}
\keyword{datagen}
