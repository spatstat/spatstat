\name{pixellate.owin}  %DontDeclareMethods
\Rdversion{1.1}
\alias{pixellate.owin}
\title{
  Convert Window to Pixel Image
}
\description{
  Convert a window to 
  a pixel image by measuring the area of intersection
  between the window and each pixel in a raster.
}
\usage{
\method{pixellate}{owin}(x, W = NULL, ...)
}
\arguments{
  \item{x}{
    Window (object of class \code{"owin"}) to be converted.
  }
  \item{W}{
    Optional. Window determining the pixel raster on which the conversion
    should occur. 
  }
  \item{\dots}{
    Optional. Extra arguments passed to \code{\link{as.mask}}
    to determine the pixel raster.
  }
}
\details{
  This is a method for the generic function \code{pixellate}.

  It converts a window \code{x} into
  a pixel image, by measuring the \emph{amount} of \code{x} that is
  inside each pixel.

  (The related function \code{\link{as.im}}
  also converts \code{x} into a pixel image,
  but records only the presence or absence of
  \code{x} in each pixel.)

  The pixel raster for the conversion is determined by the
  argument \code{W} and the extra arguments \code{\dots}.
  \itemize{
    \item 
    If \code{W} is given, and it is a binary mask (a window
    of type \code{"mask"}) then it determines the pixel raster.
    \item
    If \code{W} is given, but it is not a binary mask (it is a
    window of another type) then it will be converted to a binary mask
    using \code{as.mask(W, \dots)}.
    \item
    If \code{W} is not given, it defaults to
    \code{as.mask(as.rectangle(x), \dots)}
  }
  In the second and third cases it would be common to use the
  argument \code{dimyx} to control the number of pixels. See the Examples.

  The algorithm then computes the area
  of intersection of each pixel with the window.

  The result is a pixel image with pixel entries equal to these
  intersection areas.
}
\value{
  A pixel image (object of class \code{"im"}).
}
\seealso{
  \code{\link{pixellate.ppp}},
  \code{\link{pixellate}},
  \code{\link{as.im}}
}
\examples{
   data(letterR)
   plot(pixellate(letterR, dimyx=15))
   W <- grow.rectangle(as.rectangle(letterR), 0.2)
   plot(pixellate(letterR, W, dimyx=15))
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{manip}
