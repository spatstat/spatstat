\name{diameter.owin}  %DontDeclareMethods
\alias{diameter.owin}
\title{Diameter of a Window}
\description{
  Computes the diameter of a window.
}
\usage{
 \method{diameter}{owin}(x)
}
\arguments{
  \item{x}{
    A window whose diameter will be computed.
  }
}
\value{
  The numerical value of the diameter of the window. 
}
\details{
  This function computes the 
  diameter of a window of arbitrary shape,
  i.e. the maximum distance between any two points
  in the window.

  The argument \code{x} should be a window (an object of class
  \code{"owin"}, see \code{\link{owin.object}} for details)
  or can be given in any format acceptable to \code{\link{as.owin}()}.
  
  The function \code{diameter} is generic. This function is the
  method for the class \code{"owin"}.
}
\seealso{
  \code{\link{area.owin}},
  \code{\link{perimeter}},
  \code{\link{edges}},
  \code{\link{owin}},
  \code{\link{as.owin}}
}
\examples{
  w <- owin(c(0,1),c(0,1))
  diameter(w) 
  # returns sqrt(2)
  data(letterR)
  diameter(letterR)
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{math}
