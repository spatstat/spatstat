\name{rescue.rectangle}
\alias{rescue.rectangle}
\title{Convert Window Back To Rectangle}
\description{
  Determines whether the given window is really a rectangle
  aligned with the coordinate axes, and if so, converts it to
  a rectangle object.
}
\usage{
 rescue.rectangle(W)
}
\arguments{
  \item{W}{A window (object of class \code{"owin"}).}
}
\value{
  Another object of class \code{"owin"} representing the
  same window.
}
\details{
  This function decides whether the window \code{W} is actually a rectangle
  aligned with the coordinate axes. This will be true if
  \code{W} is 
  \itemize{
    \item a rectangle (window object of type \code{"rectangle"});
    \item a polygon (window object of type \code{"polygonal"}
    with a single polygonal boundary) that is a rectangle aligned with
    the coordinate axes;
    \item a binary mask (window object of type \code{"mask"}) in which
    all the pixel entries are \code{TRUE}.
  }
  If so, the function returns
  this rectangle, a window object of type \code{"rectangle"}.
  If not, the function returns \code{W}.
}
\seealso{
  \code{\link{as.owin}},
  \code{\link{owin.object}}
}
\examples{
  w <- owin(poly=list(x=c(0,1,1,0),y=c(0,0,1,1)))
  rw <- rescue.rectangle(w)

  w <- as.mask(unit.square())
  rw <- rescue.rectangle(w)
}
\author{Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  and Rolf Turner \email{r.turner@auckland.ac.nz}
  
}
\keyword{spatial}
\keyword{manip}
