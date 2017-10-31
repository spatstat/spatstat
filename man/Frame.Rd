\name{Frame}
\alias{Frame}
\alias{Frame<-}
\alias{Frame.default}
\alias{Frame<-.default}
\alias{Frame<-.owin}
\alias{Frame<-.ppp}
\alias{Frame<-.im}
\title{
  Extract or Change the Containing Rectangle of a Spatial Object
}
\description{
  Given a spatial object (such as a point pattern or pixel image)
  in two dimensions, these functions extract or change the
  containing rectangle inside which the object is defined.
}
\usage{
   Frame(X)

   \method{Frame}{default}(X)

   Frame(X) <- value

   \method{Frame}{owin}(X) <- value

   \method{Frame}{ppp}(X) <- value

   \method{Frame}{im}(X) <- value

   \method{Frame}{default}(X) <- value
}
\arguments{
  \item{X}{
    A spatial object such as a point pattern, line segment pattern
    or pixel image.
  }
  \item{value}{
    A rectangular window (object of class \code{"owin"}
    of type \code{"rectangle"}) to be used as the new containing
    rectangle for \code{X}.
  }
}
\details{
  The functions \code{Frame} and \code{Frame<-} are generic.

  \code{Frame(X)} extracts the rectangle inside which \code{X} is
  defined. 

  \code{Frame(X) <- R} changes the rectangle inside which \code{X} is defined
  to the new rectangle \code{R}.
}
\value{
  The result of \code{Frame} is a rectangular window (object of class
  \code{"owin"} of type \code{"rectangle"}).

  The result of \code{Frame<-} is the updated object \code{X},
  of the same class as \code{X}.
}
\author{
  \spatstatAuthors.
}
\seealso{
  \code{\link{Window}}
}
\examples{
   Frame(cells)
   X <- demopat
   Frame(X)
   Frame(X) <- owin(c(0, 11000), c(400, 8000))
}
\keyword{spatial}
\keyword{manip}
