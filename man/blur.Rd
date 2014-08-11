\name{blur}
\alias{blur}
\title{Apply Gaussian Blur to a Pixel Image}
\description{
  Applies a Gaussian blur to a pixel image.
}
\usage{
blur(x, sigma = NULL, ..., normalise=FALSE, bleed = TRUE, varcov=NULL)
}
\arguments{
  \item{x}{The pixel image. An object of class \code{"im"}.}
  \item{sigma}{
    Standard deviation of isotropic Gaussian smoothing kernel.
  }
  \item{\dots}{
    Ignored.
  }
  \item{normalise}{
    Logical flag indicating whether the output values should be divided
    by the corresponding blurred image of the window itself. See Details.
  }
  \item{bleed}{
    Logical flag indicating whether to allow blur to extend outside the
    original domain of the image. See Details.
  }
  \item{varcov}{
    Variance-covariance matrix of anisotropic Gaussian kernel.
    Incompatible with \code{sigma}.
  }
}
\details{
  This command applies a Gaussian blur to the pixel image \code{x}.

  The blurring kernel is the isotropic Gaussian kernel with standard
  deviation \code{sigma}, or the anisotropic Gaussian kernel with
  variance-covariance matrix \code{varcov}.
  The arguments \code{sigma} and \code{varcov} are incompatible.
  Also \code{sigma} may be a vector of length 2 giving the
  standard deviations of two independent Gaussian coordinates,
  thus equivalent to \code{varcov = diag(sigma^2)}.

  If the pixel values of \code{x} include some \code{NA} values
  (meaning that the image domain does not completely fill
  the rectangular frame) then these \code{NA} values are first reset to zero.

  The algorithm then computes the convolution \eqn{x \ast G}{x * G}
  of the (zero-padded) pixel
  image \eqn{x} with the specified Gaussian kernel \eqn{G}.
  
  If \code{normalise=FALSE}, then this convolution \eqn{x\ast G}{x * G}
  is returned.
  If \code{normalise=TRUE}, then the convolution \eqn{x \ast G}{x * G}
  is normalised by
  dividing it by the convolution \eqn{w \ast G}{w * G} of the image
  domain \code{w}
  with the same Gaussian kernel. Normalisation ensures that the result
  can be interpreted as a weighted average of input pixel values,
  without edge effects due to the shape of the domain.

  If \code{bleed=FALSE}, then pixel values outside the original image
  domain are set to \code{NA}. Thus the output is a pixel image with the
  same domain as the input. If \code{bleed=TRUE}, then no such
  alteration is performed, and the result is a pixel image defined
  everywhere in the rectangular frame containing the input image.
  
  Computation is performed using the Fast Fourier Transform.
}
\value{
  A pixel image with the same pixel array as the input image \code{x}.
}
\seealso{
  \code{\link{interp.im}} for interpolating a pixel image to a finer resolution,
  \code{\link{density.ppp}} for blurring a point pattern,
  \code{\link{Smooth.ppp}} for interpolating marks attached to points.
}
\examples{
   data(letterR)
   Z <- as.im(function(x,y) { 4 * x^2 + 3 * y }, letterR)
   par(mfrow=c(1,3))
   plot(Z)
   plot(letterR, add=TRUE)
   plot(blur(Z, 0.3, bleed=TRUE))
   plot(letterR, add=TRUE)
   plot(blur(Z, 0.3, bleed=FALSE))
   plot(letterR, add=TRUE)
   par(mfrow=c(1,1))
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{nonparametric}
