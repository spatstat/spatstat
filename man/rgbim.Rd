\name{rgbim}
\alias{rgbim}
\alias{hsvim}
\title{Create Colour-Valued Pixel Image}
\description{
  Creates an object of
  class \code{"im"} representing a two-dimensional pixel image
  whose pixel values are colours.
}
\usage{
  rgbim(R, G, B, maxColorValue=255)
  hsvim(H, S, V)
}
\arguments{
  \item{R,G,B}{
    Pixel images (objects of class \code{"im"}) or constants
    giving the red, green, and blue components of a colour, respectively.
  }
  \item{maxColorValue}{
    Maximum colour value for \code{R,G,B}.
  }
  \item{H,S,V}{
    Pixel images (objects of class \code{"im"}) or constants
    giving the hue, saturation, and value components of
    a colour, respectively.
  }
}
\details{
  These functions take three pixel images, with real or integer pixel values,
  and create a single pixel image whose pixel values are colours
  recognisable to \R.

  Some of the arguments may be constant numeric values, but
  at least one of the arguments must be a pixel image.
  The image arguments should be compatible (in array dimension
  and in spatial position).

  \code{rgbim} calls \code{\link{rgb}} to compute the colours,
  while \code{hsvim} calls \code{\link{hsv}}. See the help for the relevant
  function for more information about the meaning of the colour
  channels.
}
\seealso{
  \code{\link{im.object}},
  \code{\link{rgb}},
  \code{\link{hsv}}.

  See \code{\link[spatstat:colourtools]{colourtools}}
  for additional colour tools.
}
\examples{
  \testonly{
    op <- spatstat.options(npixel=32)
  }
  # create three images with values in [0,1]
  X <- setcov(owin())
  X <- eval.im(pmin(1,X))
  M <- as.owin(X)
  Y <- as.im(function(x,y){(x+1)/2}, W=M)
  Z <- as.im(function(x,y){(y+1)/2}, W=M)
  RGB <- rgbim(X, Y, Z, 1)
  HSV <- hsvim(X, Y, Z)
  plot(RGB, valuesAreColours=TRUE)
  plot(HSV, valuesAreColours=TRUE)
  \testonly{
    spatstat.options(op)
  }
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{manip}
\keyword{datagen}
 
 
