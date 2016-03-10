\name{gridcentres}
\alias{gridcentres}
\alias{gridcenters}
\title{Rectangular grid of points}
\description{
  Generates a rectangular grid of points in a window
}
\usage{
 gridcentres(window, nx, ny)
}
\arguments{
  \item{window}{A window. 
    An object of class \code{\link{owin}},
    or data in any format acceptable to \code{\link{as.owin}()}.
  }
  \item{nx}{Number of points
    in each row of the rectangular grid.
  }
  \item{ny}{Number of points
    in each column of the rectangular grid.
  }
}
\value{
  A list with two components \code{x} and \code{y}, which are numeric
  vectors giving the coordinates of the points of the
  rectangular grid.
}
\details{
  This function creates a rectangular grid of points in the window.
  
  The bounding rectangle of the \code{window} is divided into
  a regular \eqn{nx \times ny}{nx * ny} grid of rectangular tiles.
  The function returns the \eqn{x,y} coordinates of the
  centres of these tiles.

  Note that some of these grid points may lie outside the window,
  if \code{window} is not of type \code{"rectangle"}. The function
  \code{\link{inside.owin}} can be used to select those grid points
  which do lie inside the window. See the examples.

  This function is useful in creating dummy points for quadrature
  schemes (see \code{\link{quadscheme}}) and for other miscellaneous
  purposes.
}
\seealso{
  \code{\link{quad.object}},
  \code{\link{quadscheme}},
  \code{\link{inside.owin}},
  \code{\link{stratrand}}
}
\examples{
  w <- unit.square()
  xy <- gridcentres(w, 10,15)
  \dontrun{
  plot(w)
  points(xy)
  }

  bdry <- list(x=c(0.1,0.3,0.7,0.4,0.2),
               y=c(0.1,0.1,0.5,0.7,0.3))
  w <- owin(c(0,1), c(0,1), poly=bdry)
  xy <- gridcentres(w, 30, 30)
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
