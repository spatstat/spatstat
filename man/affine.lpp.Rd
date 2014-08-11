\name{affine.lpp}
\alias{affine.lpp}
\alias{shift.lpp}
\alias{rotate.lpp}
\alias{rescale.lpp}
\alias{scalardilate.lpp}
\title{Apply Geometrical Transformations to Point Pattern on a Linear Network}
\description{
  Apply geometrical transformations to a point pattern on a linear network.
}
\usage{
  \method{affine}{lpp}(X, mat=diag(c(1,1)), vec=c(0,0), \dots)

  \method{shift}{lpp}(X, \dots)

  \method{rotate}{lpp}(X, angle=pi/2, \dots)

  \method{scalardilate}{lpp}(X, f, \dots)

  \method{rescale}{lpp}(X, s)
}
\arguments{
  \item{X}{Point pattern on a linear network (object of class \code{"lpp"}).}
  \item{mat}{Matrix representing a linear transformation.}
  \item{vec}{Vector of length 2 representing a translation.}
  \item{angle}{Rotation angle in radians.}
  \item{f}{Scalar dilation factor.}
  \item{s}{
    Unit conversion factor: the new units are \code{s} times the old units.
  }
  \item{\dots}{
    Arguments passed to other methods.
  }
}
\value{
  Another point pattern on a linear network (object of class
  \code{"lpp"}) 
  representing the
  result of applying the geometrical transformation.
}
\details{
  These functions are methods for the generic functions
  \code{\link{affine}},
  \code{\link{shift}},
  \code{\link{rotate}},
  \code{\link{rescale}} and
  \code{\link{scalardilate}}
  applicable to objects of class \code{"lpp"}.

  All of these functions
  perform geometrical transformations on the object \code{X},
  except for \code{rescale}, which simply rescales the units of length.
}
\seealso{
  \code{\link{lpp}}.
  
  Generic functions
  \code{\link{affine}},
  \code{\link{shift}},
  \code{\link{rotate}},
  \code{\link{scalardilate}},
  \code{\link{rescale}}.
}
\examples{
  X <- rpoislpp(2, simplenet)
  U <- rotate(X, pi)
  stretch <- diag(c(2,3))
  Y <- affine(X, mat=stretch)
  shear <- matrix(c(1,0,0.6,1),ncol=2, nrow=2)
  Z <- affine(X, mat=shear, vec=c(0, 1))
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{math}
