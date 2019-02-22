\name{scalardilate} 
\alias{scalardilate}
\alias{scalardilate.im}
\alias{scalardilate.owin}
\alias{scalardilate.ppp}
\alias{scalardilate.psp}
\alias{scalardilate.default}
\title{Apply Scalar Dilation}
\description{
  Applies scalar dilation to a plane geometrical object,
  such as a point pattern or a window,
  relative to a specified origin.
}
\usage{
  scalardilate(X, f, \dots)

  \method{scalardilate}{im}(X, f, \dots, origin=NULL)

  \method{scalardilate}{owin}(X, f, \dots, origin=NULL)

  \method{scalardilate}{ppp}(X, f, \dots, origin=NULL)

  \method{scalardilate}{psp}(X, f, \dots, origin=NULL)

  \method{scalardilate}{default}(X, f, \dots)
}
\arguments{
  \item{X}{Any suitable dataset representing a two-dimensional
    object, such as a point pattern (object of class \code{"ppp"}),
    a window (object of class \code{"owin"}), a pixel image
    (class \code{"im"}) and so on.
  }
  \item{f}{
    Scalar dilation factor. A finite number greater than zero.
  }
  \item{\dots}{Ignored by the methods.}
  \item{origin}{
    Origin for the scalar dilation. Either a vector of 2 numbers,
    or one of the character strings
    \code{"centroid"}, \code{"midpoint"},
    \code{"left"}, \code{"right"}, \code{"top"}, \code{"bottom"},
    \code{"topleft"}, \code{"bottomleft"}, \code{"topright"} or
    \code{"bottomright"} (partially matched).
  }
}
\value{
  Another object of the same type, representing the
  result of applying the scalar dilation.
}
\details{
  This command performs scalar dilation of the object \code{X}
  by the factor \code{f} relative to the origin specified by
  \code{origin}. 

  The function \code{scalardilate} is generic, with methods for
  windows (class \code{"owin"}), point patterns (class \code{"ppp"}),
  pixel images (class \code{"im"}), line segment patterns (class \code{"psp"})
  and a default method.

  If the argument \code{origin} is not given,
  then every spatial coordinate is multiplied by the factor \code{f}.

  If \code{origin} is given, then scalar dilation is performed
  relative to the specified origin. Effectively, \code{X} is shifted
  so that \code{origin} is moved to \code{c(0,0)}, then scalar dilation
  is performed, then the result is shifted so that \code{c(0,0)} is
  moved to \code{origin}.

  This command is a special case of an affine transformation:
  see \code{\link{affine}}.
}
\seealso{
  \code{\link{affine}},
  \code{\link{shift}}
}
\examples{
plot(letterR)
plot(scalardilate(letterR, 0.7, origin="left"), col="red", add=TRUE)
}
\author{\adrian
  
  
  and \rolf
  
}
\keyword{spatial}
\keyword{math}
