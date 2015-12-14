\name{affine.im} %DontDeclareMethods
\alias{affine.im} 
\title{Apply Affine Transformation To Pixel Image}
\description{
  Applies any affine transformation of the plane (linear transformation
  plus vector shift) to a pixel image. 
}
\usage{
  \method{affine}{im}(X, mat=diag(c(1,1)), vec=c(0,0), \dots) 
}
\arguments{
  \item{X}{Pixel image (object of class \code{"im"}).}
  \item{mat}{Matrix representing a linear transformation.}
  \item{vec}{Vector of length 2 representing a translation.}
  \item{\dots}{
    Optional arguments passed to \code{\link{as.mask}}
    controlling the pixel resolution of the transformed image.
  }
}
\value{
  Another pixel image (of class \code{"im"}) representing the
  result of applying the affine transformation.
}
\details{
  The image is subjected first to the linear transformation represented by
  \code{mat} (multiplying on the left by \code{mat}),
  and then the result is translated by the vector \code{vec}.
  
  The argument \code{mat} must be a nonsingular \eqn{2 \times 2}{2 * 2}
  matrix.

  This is a method for the generic function \code{\link{affine}}.
}
\seealso{
  \code{\link{affine}},
  \code{\link{affine.ppp}},
  \code{\link{affine.psp}},
  \code{\link{affine.owin}},
  \code{\link{rotate}},
  \code{\link{shift}}
}
\examples{
  X <- setcov(owin())
  stretch <- diag(c(2,3))
  Y <- affine(X, mat=stretch)
  shear <- matrix(c(1,0,0.6,1),ncol=2, nrow=2)
  Z <- affine(X, mat=shear)
}
\author{Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  and Rolf Turner \email{r.turner@auckland.ac.nz}
  
}
\keyword{spatial}
\keyword{math}
