\name{interp.im}
\alias{interp.im}
\title{Interpolate a Pixel Image}
\description{
  Interpolates the values of a pixel image at any
  desired location in the frame.
}
\usage{
interp.im(Z, x, y)
}
\arguments{
  \item{Z}{
    Pixel image (object of class \code{"im"})
    with numeric or integer values.
  }
  \item{x,y}{
    Vectors of Cartesian coordinates.
  }
}
\details{
  A value at each location \code{(x[i],y[i])} will be
  interpolated using the pixel values of \code{Z} at the four
  surrounding pixel centres, by simple bilinear interpolation.

  At the boundary (where \code{(x[i],y[i])} is not surrounded by
  four pixel centres) the value at the nearest pixel
  is taken.
}
\value{
  Vector of interpolated values, with \code{NA} for points that lie
  outside the domain of the image.
}
\examples{
   opa <- par(mfrow=c(1,2))
   # coarse image
   V <- as.im(function(x,y) { x^2 + y }, owin(), dimyx=10)
   plot(V, main="coarse image", col=terrain.colors(256))

   # lookup value at location (0.5,0.5)
   V[list(x=0.5,y=0.5)]
   # interpolated value at location (0.5,0.5)
   interp.im(V, 0.5, 0.5)
   # true value is 0.75

   # how to obtain an interpolated image at a desired resolution
   U <- as.im(interp.im, W=owin(), Z=V, dimyx=256)
   plot(U, main="interpolated image", col=terrain.colors(256))
   par(opa)
}

\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{manip}
