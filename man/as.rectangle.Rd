\name{as.rectangle}
\alias{as.rectangle}
\title{Window Frame}
\description{
  Extract the window frame of a window
  or other spatial dataset
}
\usage{
 as.rectangle(w, \dots)
}
\arguments{
  \item{w}{A window, or a dataset that has a window.
    Either a window (object of class \code{"owin"}),
    a pixel image (object of class \code{"im"})
    or other data determining such a window.
  }
  \item{\dots}{
    Optional. Auxiliary data to help determine the window. If
    \code{w} does not belong to a recognised class,
    the arguments \code{w} and \code{\dots} 
    are passed to \code{\link{as.owin}} to determine the window.
  }
}
\value{
  A window (object of class \code{"owin"})
  of type \code{"rectangle"} representing a rectangle.
}
\details{
  This function is the quickest way to determine a bounding rectangle
  for a spatial dataset.

  If \code{w} is a window, the function
  just extracts the outer bounding rectangle
  of \code{w} as given by its elements \code{xrange,yrange}.

  The function can also be applied to any spatial dataset that has a window:
  for example, a point pattern (object of class \code{"ppp"}) or
  a line segment pattern (object of class \code{"psp"}).
  The bounding rectangle of the window of the dataset is extracted.

  Use the function \code{\link{boundingbox}} to compute the \emph{smallest}
  bounding rectangle of a dataset.
}
\seealso{
  \code{\link{owin}},
  \code{\link{as.owin}},
  \code{\link{boundingbox}}
}
\examples{
  w <- owin(c(0,10),c(0,10), poly=list(x=c(1,2,3,2,1), y=c(2,3,4,6,7)))
  r <- as.rectangle(w)
  # returns a 10 x 10 rectangle

  data(lansing)
  as.rectangle(lansing)

  data(copper)
  as.rectangle(copper$SouthLines)
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{manip}
