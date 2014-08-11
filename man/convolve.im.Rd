\name{convolve.im}
\alias{convolve.im}
\title{Convolution of Pixel Images}
\description{
  Computes the convolution of two pixel images.
}
\usage{
 convolve.im(X, Y=X, \dots, reflectX=FALSE, reflectY=FALSE)
}
\arguments{
  \item{X}{
    A pixel image (object of class \code{"im"}.
  }
  \item{Y}{
    Optional. Another pixel image.
  }
  \item{\dots}{Ignored.}
  \item{reflectX,reflectY}{
    Logical values specifying whether the images \code{X} and \code{Y}
    (respectively) should be reflected in the origin before computing
    the convolution.
  }
}
\value{
  A pixel image (an object of class \code{"im"}) representing the
  convolution of \code{X} and \code{Y}.
}
\details{
  The \emph{convolution} of two pixel images \eqn{X} and \eqn{Y} in the plane
  is the function \eqn{C(v)} defined for each vector \eqn{v} as
  \deqn{
    C(v) = \int X(u)Y(v-u)\, {\rm d}u
  }{
    C(v) = integral of X(u) * Y(v-u) du
  }
  where the integral is
  over all spatial locations \eqn{u}, and where \eqn{X(u)} and
  \eqn{Y(u)} denote the pixel values of \eqn{X} and \eqn{Y} respectively
  at location \eqn{u}.
  
  This command computes a discretised approximation to
  the convolution, using the Fast Fourier Transform.
  The return value is
  another pixel image (object of class \code{"im"}) whose greyscale values
  are values of the convolution.

  If \code{reflectX = TRUE} then the pixel image \code{X} is reflected
  in the origin (see \code{\link{reflect}})
  before the convolution is computed, so that
  \code{convolve.im(X,Y,reflectX=TRUE)} is mathematically
  equivalent to \code{convolve.im(reflect(X), Y)}. (These two commands
  are not exactly equivalent, because the reflection is performed
  in the Fourier domain in the first command, and reflection is
  performed in the spatial domain in the second command).

  Similarly if \code{reflectY = TRUE} then the pixel image \code{Y} is reflected
  in the origin before the convolution is computed, so that
  \code{convolve.im(X,Y,reflectY=TRUE)} is mathematically
  equivalent to \code{convolve.im(X, reflect(Y))}.
}
\seealso{
  \code{\link{imcov}},
  \code{\link{reflect}}
}
\examples{
  X <- as.im(letterR)
  Y <- as.im(square(1))
  plot(convolve.im(X, Y))
  plot(convolve.im(X, Y, reflectX=TRUE))
  plot(convolve.im(X))
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{math}
