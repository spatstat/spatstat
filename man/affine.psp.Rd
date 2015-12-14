\name{affine.psp}  %DontDeclareMethods
\alias{affine.psp}
\title{Apply Affine Transformation To Line Segment Pattern}
\description{
  Applies any affine transformation of the plane (linear transformation
  plus vector shift) to a line segment pattern. 
}
\usage{
 \method{affine}{psp}(X, mat=diag(c(1,1)), vec=c(0,0), \dots)
}
\arguments{
  \item{X}{Line Segment pattern (object of class \code{"psp"}).}
  \item{mat}{Matrix representing a linear transformation.}
  \item{vec}{Vector of length 2 representing a translation.}
  \item{\dots}{Arguments passed to \code{\link{affine.owin}} affecting
    the handling of the observation window, if it is a binary pixel
    mask.
  }
}
\value{
  Another line segment pattern (of class \code{"psp"}) representing the
  result of applying the affine transformation.
}
\details{
  The line segment pattern, and its window, are subjected first to the
  linear transformation represented by
  \code{mat} (multiplying on the left by \code{mat}),
  and are then translated by the vector \code{vec}.
  
  The argument \code{mat} must be a nonsingular \eqn{2 \times 2}{2 * 2}
  matrix.

  This is a method for the generic function \code{\link{affine}}.
}
\seealso{
  \code{\link{affine}},
  \code{\link{affine.owin}},
  \code{\link{affine.ppp}},
  \code{\link{affine.im}},
  \code{\link{flipxy}},
  \code{\link{rotate}},
  \code{\link{shift}}
}
\examples{
  oldpar <- par(mfrow=c(2,1))
  X <- psp(runif(10), runif(10), runif(10), runif(10), window=owin())
  plot(X, main="original")
  # shear transformation
  Y <- affine(X, matrix(c(1,0,0.6,1),ncol=2))
  plot(Y, main="transformed")
  par(oldpar)
  # 
  # rescale y coordinates by factor 0.2
  affine(X, diag(c(1,0.2)))
}
\author{Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  and Rolf Turner \email{r.turner@auckland.ac.nz}
  
}
\keyword{spatial}
\keyword{math}
