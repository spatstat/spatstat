\name{nearest.raster.point}
\alias{nearest.raster.point}
\title{Find Pixel Nearest to a Given Point}
\description{
  Given cartesian coordinates, find the nearest pixel.
}
\usage{
 nearest.raster.point(x,y,w, indices=TRUE)
}
\arguments{
  \item{x}{Numeric vector of \eqn{x} coordinates of any points}
  \item{y}{Numeric vector of \eqn{y} coordinates of any points}
  \item{w}{An image (object of class \code{"im"})
    or a binary mask window (an object of class \code{"owin"}
    of type \code{"mask"}).
  }
  \item{indices}{Logical flag indicating whether to return the
    row and column indices, or the actual \eqn{x,y} coordinates.
    }
}
\value{
  If \code{indices=TRUE}, a 
  list containing two vectors \code{rr} and \code{cc}
  giving row and column positions (in the image matrix).
  If \code{indices=FALSE}, a list containing
  vectors \code{x} and \code{y} giving actual coordinates
  of the pixels.
}
\details{
  The argument \code{w} should be either a pixel image
  (object of class \code{"im"}) or a window (an object of class
  \code{"owin"}, see \code{\link{owin.object}} for details)
  of type \code{"mask"}. 

  The arguments \code{x} and \code{y} should be numeric vectors
  of equal length. They are interpreted as the coordinates of
  points in space. For each point \code{(x[i], y[i])}, the function
  finds the nearest pixel in the grid of pixels for \code{w}.

  If \code{indices=TRUE},
  this function returns a list containing two vectors \code{rr} and
  \code{cc} giving row and column positions (in the image matrix).
  For the location \code{(x[i],y[i])} the nearest
  pixel is at row \code{rr[i]} and column \code{cc[i]} of
  the image.

  If \code{indices=FALSE}, the function returns a list containing
  two vectors \code{x} and \code{y} giving the actual coordinates
  of the pixels.
}
\seealso{
  \code{\link{owin.object}},
  \code{\link{as.mask}}
}
\examples{
  w <- owin(c(0,1), c(0,1), mask=matrix(TRUE, 100,100))  # 100 x 100 grid
  nearest.raster.point(0.5, 0.3, w)
  nearest.raster.point(0.5, 0.3, w, indices=FALSE)
}
\author{Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  and Rolf Turner \email{r.turner@auckland.ac.nz}
  
}
\keyword{spatial}
\keyword{manip}
