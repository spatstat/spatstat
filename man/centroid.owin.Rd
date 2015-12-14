\name{centroid.owin}
\alias{centroid.owin}
\title{Centroid of a window}
\description{
  Computes the centroid (centre of mass) of a window 
}
\usage{
 centroid.owin(w, as.ppp = FALSE)
}
\arguments{
  \item{w}{A window}
  \item{as.ppp}{Logical flag indicating whether to return the centroid
    as a point pattern (\code{ppp} object)}
}
\value{
  Either a list with components \code{x, y}, or a point pattern (of class
  \code{ppp}) consisting of a single point, giving the coordinates of the
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
  the window, unless the window is convex.
  If  \code{as.ppp=TRUE} and
  the centroid of \code{w} lies outside \code{w},
  then the window of the returned point pattern
  will be a rectangle containing the
  original window (using \code{\link{as.rectangle}}. 
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
  w <- Window(demopat)
  # an irregular window
  cent <- centroid.owin(w, as.ppp = TRUE)
  \dontrun{
  plot(cent)
  # plot the window and its centroid
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
\author{Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  ,
  Rolf Turner \email{r.turner@auckland.ac.nz}
  
  and Ege Rubak \email{rubak@math.aau.dk}
}
\keyword{spatial}
\keyword{math}




