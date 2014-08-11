\name{affine.ppp} %DontDeclareMethods
\alias{affine.ppp} 
\title{Apply Affine Transformation To Point Pattern}
\description{
  Applies any affine transformation of the plane (linear transformation
  plus vector shift) to a point pattern. 
}
\usage{
 \method{affine}{ppp}(X, mat=diag(c(1,1)), vec=c(0,0), \dots)
}
\arguments{
  \item{X}{Point pattern (object of class \code{"ppp"}).}
  \item{mat}{Matrix representing a linear transformation.}
  \item{vec}{Vector of length 2 representing a translation.}
  \item{\dots}{Arguments passed to \code{\link{affine.owin}} affecting
    the handling of the observation window, if it is a binary pixel
    mask.
  }
}
\value{
  Another point pattern (of class \code{"ppp"}) representing the
  result of applying the affine transformation.
}
\details{
  The point pattern, and its window, are subjected first to the
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
  \code{\link{affine.psp}},
  \code{\link{affine.im}},
  \code{\link{flipxy}},
  \code{\link{rotate}},
  \code{\link{shift}}
}
\examples{
  data(cells)
  # shear transformation
  X <- affine(cells, matrix(c(1,0,0.6,1),ncol=2))
  \dontrun{
  plot(X)
  # rescale y coordinates by factor 1.3
  plot(affine(cells, diag(c(1,1.3))))
  }
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{math}
