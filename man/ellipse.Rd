\name{ellipse}
\alias{ellipse}
\title{
  Elliptical Window.
}
\description{
  Create an elliptical window.
}
\usage{
   ellipse(a, b, centre=c(0,0), phi=0, \dots, mask=FALSE, npoly = 128)
}
\arguments{
  \item{a,b}{
    The half-lengths of the axes of the ellipse.
  }
  \item{centre}{
    The centre of the ellipse.
  }
  \item{phi}{
    The (anti-clockwise) angle through which the ellipse should be
    rotated (about its centre) starting from an orientation in which
    the axis of half-length \code{a} is horizontal.
  }
  \item{mask}{
    Logical value controlling the type of approximation
    to a perfect ellipse. See Details.
  }
  \item{\dots}{
    Arguments passed to \code{\link{as.mask}} to determine the
    pixel resolution, if \code{mask} is \code{TRUE}.
  }
  \item{npoly}{
    The number of edges in the polygonal approximation to the ellipse.
  }
}
\details{
  This command creates a window object
  representing an ellipse with the given centre and axes.

  By default, the ellipse is
  approximated by a polygon with \code{npoly} edges.

  If \code{mask=TRUE}, then the ellipse is approximated by a binary pixel
  mask. The resolution of the mask is controlled by
  the arguments \code{\dots} which are passed to \code{\link{as.mask}}.

  The arguments \code{a} and \code{b} must be single positive numbers.
  The argument \code{centre} specifies the ellipse centre: it can be either 
  a numeric vector of length 2 giving the coordinates,
  or a \code{list(x,y)} giving the coordinates of exactly one point, or a
  point pattern (object of class \code{"ppp"}) containing exactly one point.
}
\value{
  An object of class \code{owin} (either of type \dQuote{polygonal}
  or of type \dQuote{mask}) specifying an elliptical window.
}
\author{\adrian
  
  
  and \rolf
  
}
\seealso{
  \code{\link{disc}},
  \code{\link{owin.object}},
  \code{\link{owin}},
  \code{\link{as.mask}}
}
\examples{
  W <- ellipse(a=5,b=2,centre=c(5,1),phi=pi/6)
  plot(W,lwd=2,border="red")
  WM <- ellipse(a=5,b=2,centre=c(5,1),phi=pi/6,mask=TRUE,dimyx=512)
  plot(WM,add=TRUE,box=FALSE)
}
\keyword{spatial}
\keyword{datagen}
