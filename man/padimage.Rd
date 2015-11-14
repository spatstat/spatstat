\name{padimage}
\alias{padimage}
\title{
  Pad the Border of a Pixel Image
}
\description{
  Fills the border of a pixel image with a given value or values,
  or extends a pixel image to fill a larger window.
}
\usage{
padimage(X, value=NA, n=1, W=NULL)
}
\arguments{
  \item{X}{
    Pixel image (object of class \code{"im"}).
  }
  \item{value}{
    Single value to be placed around the border of \code{X}.    
  }
  \item{n}{
    Width of border, in pixels.
    See Details.
  }
  \item{W}{
    Window for the resulting image.
    Incompatible with \code{n}.
  }
}
\details{
  The image \code{X} will be expanded by a margin of \code{n}
  pixels, or extended to fill the window \code{W},
  with new pixel values set to \code{value}.

  The argument \code{value} should be a single value (a vector of length
  1), normally a value of the same type as the pixel values of \code{X}.
  It may be \code{NA}. Alternatively if \code{X} is a
  factor-valued image, \code{value} can be one of the levels of
  \code{X}.

  If \code{n} is given, it may be a single number, specifying the
  width of the border in pixels. Alternatively it may be a vector
  of length 2 or 4. It will be replicated to length 4, and these numbers
  will be interpreted as the border widths for the
  (left, right, top, bottom) margins respectively.

  Alternatively if \code{W} is given, the image will be extended to
  the window \code{W}. 
}
\value{
  Another object of class \code{"im"}, of the same type as \code{X}.
}
\author{\adrian
  
  
  \rolf
  
  and \ege
  
}
\seealso{
  \code{\link{im}}
}
\examples{
 Z <- setcov(owin())
 plot(padimage(Z, 1, 10))
}
\keyword{spatial}
\keyword{manip}

