\name{trim.rectangle}
\alias{trim.rectangle}
\title{Cut margins from rectangle}
\description{
  Trims a margin from a rectangle.
}
\usage{
 trim.rectangle(W, xmargin=0, ymargin=xmargin)
}
\arguments{
  \item{W}{
    A window (object of class \code{"owin"}).
    Must be of type \code{"rectangle"}.
  }
  \item{xmargin}{Width of horizontal margin to be trimmed.
    A single nonnegative number, or a vector of length 2
    indicating margins of unequal width at left and right.
  }
  \item{ymargin}{Height of vertical margin to be trimmed.
    A single nonnegative number, or a vector of length 2
    indicating margins of unequal width at bottom and top.
  }
}
\value{
  Another object of class \code{"owin"} representing the
  window after margins are trimmed.
}
\details{
  This is a simple convenience function to trim off a
  margin of specified width and height from each side of a
  rectangular window. Unequal margins can also be trimmed.
}
\seealso{
  \code{\link{owin.object}}
}
\examples{
  w <- square(10)
  # trim a margin of width 1 from all four sides
  square9 <- trim.rectangle(w, 1)

  # trim margin of width 3 from the right side
  # and margin of height 4 from top edge.
  v <- trim.rectangle(w, c(0,3), c(0,4))
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{manip}
