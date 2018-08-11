\name{boundingbox}
\alias{boundingbox}
\alias{boundingbox.default}
\alias{boundingbox.im}
\alias{boundingbox.owin}
\alias{boundingbox.ppp}
\alias{boundingbox.psp}
\alias{boundingbox.lpp}
\alias{boundingbox.linnet}
\alias{boundingbox.solist}
\title{
  Bounding Box of a Window, Image, or Point Pattern
}
\description{
  Find the smallest rectangle containing a given window(s),
  image(s) or point pattern(s).
}
\usage{
boundingbox(\dots)

\method{boundingbox}{default}(\dots)

\method{boundingbox}{im}(\dots)

\method{boundingbox}{owin}(\dots)

\method{boundingbox}{ppp}(\dots)

\method{boundingbox}{psp}(\dots)

\method{boundingbox}{lpp}(\dots)

\method{boundingbox}{linnet}(\dots)

\method{boundingbox}{solist}(\dots)
}
\arguments{
  \item{\dots}{One or more windows (objects of class \code{"owin"}),
    pixel images (objects of class \code{"im"}) or
    point patterns (objects of class \code{"ppp"} or \code{"lpp"})
    or line segment patterns (objects of class \code{"psp"})
    or linear networks (objects of class \code{"linnet"})
    or any combination of such objects.
    Alternatively, the argument may be a list of such objects,
    of class \code{"solist"}.
  }
}
\details{
  This function finds the smallest rectangle (with sides parallel to
  the coordinate axes) that contains all the given objects.

  For a window (object of class \code{"owin"}), the bounding box
  is the smallest rectangle that contains all the vertices of the
  window (this is generally smaller than the enclosing frame,
  which is returned by \code{\link{as.rectangle}}).

  For a point pattern (object of class \code{"ppp"} or \code{"lpp"}),
  the bounding box
  is the smallest rectangle that contains all the points of the pattern.
  This is usually smaller than the bounding box of the window of the
  point pattern.

  For a line segment pattern (object of class \code{"psp"})
  or a linear network (object of class \code{"linnet"}), the
  bounding box is the smallest rectangle that contains all endpoints
  of line segments.
  
  For a pixel image (object of class \code{"im"}), the image will
  be converted to a window using \code{\link{as.owin}},
  and the bounding box of this window is obtained.
  
  If the argument is a list of several objects, then
  this function finds the smallest rectangle that contains
  all the bounding boxes of the objects.
}
\value{
  \code{\link{owin}},
  \code{\link{as.owin}},
  \code{\link{as.rectangle}}
}
\examples{
  w <- owin(c(0,10),c(0,10), poly=list(x=c(1,2,3,2,1), y=c(2,3,4,6,7)))
  r <- boundingbox(w)
  # returns rectangle [1,3] x [2,7]

  w2 <- unit.square()
  r <- boundingbox(w, w2)
  # returns rectangle [0,3] x [0,7]
}
\author{
  \spatstatAuthors.
}
\keyword{spatial}
\keyword{utilities}
