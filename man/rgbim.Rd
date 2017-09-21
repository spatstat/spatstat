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
  rgbim(R, G, B, A, maxColorValue=255, autoscale=FALSE)
  hsvim(H, S, V, A, autoscale=FALSE)
}
\arguments{
  \item{R,G,B}{
    Pixel images (objects of class \code{"im"}) or constants
    giving the red, green, and blue components of a colour, respectively.
  }
  \item{A}{
    Optional. Pixel image or constant value
    giving the alpha (transparency) component of a colour.
  }
  \item{maxColorValue}{
    Maximum colour channel value for \code{R,G,B,A}.
  }
  \item{H,S,V}{
    Pixel images (objects of class \code{"im"}) or constants
    giving the hue, saturation, and value components of
    a colour, respectively.
  }
  \item{autoscale}{
    Logical. If \code{TRUE}, input values are automatically rescaled
    to fit the permitted range. RGB values are scaled to lie between
    0 and \code{maxColorValue}. HSV values are scaled to lie between 0
    and 1.
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
  M <- Window(X)
  Y <- as.im(function(x,y){(x+1)/2}, W=M)
  Z <- as.im(function(x,y){(y+1)/2}, W=M)
  RGB <- rgbim(X, Y, Z, maxColorValue=1)
  HSV <- hsvim(X, Y, Z)
  plot(RGB, valuesAreColours=TRUE)
  plot(HSV, valuesAreColours=TRUE)
  \testonly{
    spatstat.options(op)
  }
}
\author{\adrian
  
  
  and \rolf
  
}
\keyword{spatial}
\keyword{manip}
\keyword{datagen}
 
 
