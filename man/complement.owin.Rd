\name{complement.owin}
\alias{complement.owin}
\title{Take Complement of a Window}
\description{
  Take the set complement of a window, within its enclosing rectangle
  or in a larger rectangle.
}
\usage{
 complement.owin(w, frame=as.rectangle(w))
}
\arguments{
  \item{w}{
    an object of class \code{"owin"} describing
    a window of observation for a point pattern.
  }
  \item{frame}{
    Optional. The enclosing rectangle,
    with respect to which the set complement is taken.
  }
}
\value{
  Another object of class \code{"owin"} 
  representing the complement of the window, i.e. the inside
  of the window becomes the outside.
}
\details{
  This yields a window object (of class \code{"owin"},
  see \code{\link{owin.object}}) representing the set complement
  of \code{w} with respect to the rectangle \code{frame}.

  By default, \code{frame} is the enclosing box of \code{w}
  (originally specified by the arguments \code{xrange} and \code{yrange}
  given to \code{\link{owin}} when \code{w} was created).
  If \code{frame} is specified, it must be a rectangle (an object of
  class \code{"owin"} whose type is \code{"rectangle"}) and it must be
  larger than the enclosing box of \code{w}. This rectangle becomes the
  enclosing box for the resulting window.
  
  If \code{w} is a rectangle, then \code{frame} must be specified.
  Otherwise an error will occur (since the complement of \code{w} in
  itself is empty).

  For rectangular and polygonal windows, the complement is computed by
  reversing the sign of each boundary polygon, while for binary masks it is
  computed by negating the pixel values.
}
\seealso{
  \code{\link{owin}},
  \code{\link{owin.object}}
}
\examples{
   # rectangular
   a <- owin(c(0,1),c(0,1))
   b <- owin(c(-1,2),c(-1,2))
   bmina <- complement.owin(a, frame=b)
   # polygonal
   data(demopat)
   w <- demopat$window
   outside <- complement.owin(w)
   # mask
   w <- as.mask(demopat$window)
   outside <- complement.owin(w)
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{math}
