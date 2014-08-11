\name{corners}
\alias{corners}
\title{Corners of a rectangle}
\description{
  Returns the four corners of a rectangle 
}
\usage{
 corners(window)
}
\arguments{
  \item{window}{A window. 
    An object of class \code{\link{owin}},
    or data in any format acceptable to \code{\link{as.owin}()}.
  }
}
\value{
  A list with two components \code{x} and \code{y}, which are numeric
  vectors of length 4 giving the coordinates of the four corner points
  of the (bounding rectangle of the) window.
}
\details{
  This trivial function is occasionally convenient.
  If \code{window} is of type \code{"rectangle"} this returns the
  four corners of the window itself; 
  otherwise, it returns the corners of the bounding rectangle
  of the window.
}
\seealso{
  \code{\link{quad.object}},
  \code{\link{quadscheme}}
}
\examples{
  w <- unit.square()
  corners(w)
       # returns list(x=c(0,1,0,1),y=c(0,0,1,1))
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{utilities}
