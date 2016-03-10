\name{affine}
\alias{affine}
\title{Apply Affine Transformation}
\description{
  Applies any affine transformation of the plane (linear transformation
  plus vector shift) to a plane geometrical object,
  such as a point pattern or a window. 
}
\usage{
  affine(X, \dots)
}
\arguments{
  \item{X}{Any suitable dataset representing a two-dimensional
    object, such as a point pattern (object of class \code{"ppp"}),
    a line segment pattern (object of class \code{"psp"}),
    a window (object of class \code{"owin"}) or a pixel image
    (object of class \code{"im"}).
  }
  \item{\dots}{Arguments determining the affine transformation.}
}
\value{
  Another object of the same type, representing the
  result of applying the affine transformation.
}
\details{
  This is generic. Methods are provided for
  point patterns (\code{\link{affine.ppp}})
  and windows (\code{\link{affine.owin}}).
}
\seealso{
  \code{\link{affine.ppp}},
  \code{\link{affine.psp}},
  \code{\link{affine.owin}},
  \code{\link{affine.im}},
  \code{\link{flipxy}},
  \code{\link{reflect}},
  \code{\link{rotate}},
  \code{\link{shift}}
}
\author{Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  and Rolf Turner \email{r.turner@auckland.ac.nz}
  
}
\keyword{spatial}
\keyword{math}
