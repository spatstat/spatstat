\name{affine.linnet} %DontDeclareMethods
\alias{affine.linnet} 
\alias{shift.linnet}
\alias{rotate.linnet}
\alias{rescale.linnet}
\alias{scalardilate.linnet}
\title{Apply Geometrical Transformations to a Linear Network}
\description{
  Apply geometrical transformations
  to a linear network.
}
\usage{
  \method{affine}{linnet}(X, mat=diag(c(1,1)), vec=c(0,0), \dots)

  \method{shift}{linnet}(X, vec=c(0,0), \dots, origin=NULL)

  \method{rotate}{linnet}(X, angle=pi/2, \dots, centre=NULL)

  \method{scalardilate}{linnet}(X, f, \dots)

  \method{rescale}{linnet}(X, s, unitname)
}
\arguments{
  \item{X}{Linear network (object of class \code{"linnet"}).}
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
  \item{origin}{
    Character string determining a location
    that will be shifted to the origin. Options are
    \code{"centroid"}, \code{"midpoint"} and \code{"bottomleft"}.
    Partially matched. 
  }
  \item{centre}{
    Centre of rotation.
    Either a vector of length 2, or a character string
    (partially matched to \code{"centroid"}, \code{"midpoint"}
    or \code{"bottomleft"}).
    The default is the coordinate origin \code{c(0,0)}.
  }
  \item{unitname}{
    Optional. New name for the unit of length.
    A value acceptable to the function \code{\link{unitname<-}}
  }
}
\value{
  Another linear network (of class \code{"linnet"}) representing the
  result of applying the geometrical transformation.
}
\details{
  These functions are methods for the generic functions
  \code{\link{affine}},
  \code{\link{shift}},
  \code{\link{rotate}},
  \code{\link{rescale}} and
  \code{\link{scalardilate}}
  applicable to objects of class \code{"linnet"}.

  All of these functions
  perform geometrical transformations on the object \code{X},
  except for \code{rescale}, which simply rescales the units of length.
}
\seealso{
  \code{\link{linnet}} and \code{\link{as.linnet}}.
  
  Generic functions
  \code{\link{affine}},
  \code{\link{shift}},
  \code{\link{rotate}},
  \code{\link{scalardilate}},
  \code{\link{rescale}}.
}
\examples{
  U <- rotate(simplenet, pi)
  stretch <- diag(c(2,3))
  Y <- affine(simplenet, mat=stretch)
  shear <- matrix(c(1,0,0.6,1),ncol=2, nrow=2)
  Z <- affine(simplenet, mat=shear, vec=c(0, 1))
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{math}
