\name{as.mask}
\alias{as.mask}
\title{Pixel Image Approximation of a Window}
\description{
  Obtain a discrete (pixel image) approximation of a given window
}
\usage{
 as.mask(w, eps=NULL, dimyx=NULL, xy=NULL)
}
\arguments{
  \item{w}{A window (object of class \code{"owin"}) or data acceptable
    to \code{\link{as.owin}}.}
  \item{eps}{(optional) width and height of pixels.}
  \item{dimyx}{(optional) pixel array dimensions}
  \item{xy}{(optional) data containing pixel coordinates}
}
\value{
  A window (object of class \code{"owin"})
  of type \code{"mask"} representing a binary pixel image.
}
\details{
  This function generates a rectangular grid of locations in the plane,
  tests whether each of these locations lies inside the
  window \code{w}, and stores the results as a binary pixel image
  or `mask' (an object of class \code{"owin"}, see \code{\link{owin.object}}).

  The most common use of this function is to approximate the shape
  of another window \code{w} by a binary pixel image. In this case,
  we will usually want to have a very fine grid of pixels.

  This function can also be used to generate a coarsely-spaced grid of
  locations inside a window, for purposes such as subsampling
  and prediction.

  The grid spacing and location are controlled by the
  arguments \code{eps}, \code{dimyx} and \code{xy},
  which are mutually incompatible.

  If \code{eps} is given, then it determines the grid spacing.
  If \code{eps} is a single number,
  then the grid spacing will be approximately \code{eps}
  in both the \eqn{x} and \eqn{y} directions. If \code{eps} is a
  vector of length 2, then the grid spacing will be approximately
  \code{eps[1]} in the \eqn{x} direction and 
  \code{eps[2]} in the \eqn{y} direction.

  If \code{dimyx} is given, then the pixel grid will be an
  \eqn{m \times n}{m x n} rectangular grid
  where \eqn{m, n} are given by \code{dimyx[2]}, \code{dimyx[1]}
  respectively. \bold{Warning:} \code{dimyx[1]} is the number of
  pixels in the \eqn{y} direction, and \code{dimyx[2]} is the number
  in the \eqn{x} direction. 

  If \code{xy} is given, then this should be some kind of
  data specifing the coordinates of a pixel grid. It may be
  \itemize{
    \item 
    a list or structure containing elements \code{x} and \code{y}
    which are numeric vectors of equal length. These will be taken as
    \eqn{x} and \code{y} coordinates of the margins
    of the grid. The pixel coordinates will be generated
    from these two vectors.
    \item
    a pixel image (object of class \code{"im"}).
    \item
    a window (object of class \code{"owin"}) which is
    of type \code{"mask"} so that it contains pixel coordinates.
  }
  If \code{xy} is given, \code{w} may be omitted.
  
  If neither \code{eps} nor \code{dimyx} nor \code{xy} is given,
  the pixel raster dimensions are obtained from
  \code{\link{spatstat.options}("npixel")}.

  There is no inverse of this function. However, the function
  \code{\link{as.polygonal}} will compute a polygonal approximation
  of a binary mask.
}
\seealso{
  \code{\link{owin.object}},
  \code{\link{as.rectangle}},
  \code{\link{as.polygonal}},
  \code{\link{spatstat.options}}
}
\examples{
  w <- owin(c(0,10),c(0,10), poly=list(x=c(1,2,3,2,1), y=c(2,3,4,6,7)))
  \dontrun{plot(w)}
  m <- as.mask(w)
  \dontrun{plot(m)}
  x <- 1:9
  y <- seq(0.25, 9.75, by=0.5)
  m <- as.mask(w, xy=list(x=x, y=y))
}
\author{Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  and Rolf Turner \email{r.turner@auckland.ac.nz}
  
}
\keyword{spatial}
\keyword{manip}
