\name{pixelcentres}
\alias{pixelcentres}
\title{
  Extract Pixel Centres as Point Pattern
}
\description{
  Given a pixel image or binary mask window,
  extract the centres of all pixels and return them as a point pattern.
}
\usage{
pixelcentres(X, W = NULL, ...)
}
\arguments{
  \item{X}{
    Pixel image (object of class \code{"im"})
    or window (object of class \code{"owin"}).
  }
  \item{W}{
    Optional window to contain the resulting point pattern.
  }
  \item{\dots}{
    Optional arguments defining the pixel resolution.
  }
}
\details{
  If the argument \code{X} is a pixel image, the result is a point
  pattern, consisting of the centre of every pixel whose pixel value
  is not \code{NA}. 

  If \code{X} is a window which is a binary mask, the result is a point
  pattern consisting of the centre of every pixel inside the window
  (i.e. every pixel for which the mask value is \code{TRUE}).

  Otherwise, \code{X} is first converted to a window, then
  converted to a mask using \code{\link{as.mask}}, then handled as
  above.
}
\value{
  A point pattern (object of class \code{"ppp"}).
}
\seealso{
  \code{\link{raster.xy}}
}
\examples{
  pixelcentres(letterR, dimyx=5)
}
\author{\adrian
  
  ,
  \rolf
  
  and \ege
}
\keyword{spatial}
\keyword{manip}
