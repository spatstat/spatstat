\name{centroid.owin}
\alias{centroid.owin}
\title{Centroid of a window}
\description{
  Computes the centroid (centre of mass) of a window 
}
\usage{
 centroid.owin(w)
}
\arguments{
  \item{w}{A window}
}
\value{
  A list with components \code{x, y} giving the coordinates of the
  centroid of the window \code{w}.
}
\details{
  The centroid of the window \code{w} is computed.
  The centroid (``centre of mass'') 
  is the point whose \eqn{x} and \eqn{y} coordinates 
  are the mean values of the \eqn{x} and \eqn{y} coordinates
  of all points in the window.

  The argument \code{w} should be a window (an object of class
  \code{"owin"}, see \code{\link{owin.object}} for details)
  or can be given in any format acceptable to \code{\link{as.owin}()}.

  The calculation uses an exact analytic formula for the case
  of polygonal windows.

  Note that the centroid of a window is not necessarily inside 
  the window. If the window is convex then it does contain its centroid.
}
\seealso{
  \code{\link{owin}},
  \code{\link{as.owin}}
}
\examples{
  w <- owin(c(0,1),c(0,1))
  centroid.owin(w) 
  # returns 0.5, 0.5

  data(demopat)
  w <- demopat$window
  # an irregular window
  \dontrun{
  plot(w)
  # plot the window
  points(centroid.owin(w))
  # mark its centroid 
  }

  wapprox <- as.mask(w)
  # pixel approximation of window
  \dontrun{
  points(centroid.owin(wapprox))
  # should be indistinguishable 
  }
  \testonly{
  centroid.owin(w)
  centroid.owin(wapprox)
  }	
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{math}




