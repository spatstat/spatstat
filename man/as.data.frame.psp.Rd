\name{as.data.frame.psp}
\alias{as.data.frame.psp}
\title{Coerce Line Segment Pattern to a Data Frame}
\description{
  Extracts the coordinates of the endpoints in a line segment pattern,
  and their marks if any, and returns them in a data frame.
}
\usage{
\method{as.data.frame}{psp}(x, row.names = NULL, ...)
}
\arguments{
  \item{x}{Line segment pattern (object of class \code{"psp"}).}
  \item{row.names}{Optional character vector of row names.}
  \item{\dots}{Ignored.}
}
\details{
  This is a method for the generic function \code{\link{as.data.frame}}
  for the class \code{"psp"} of line segment patterns.

  It extracts the coordinates of the endpoints of the line segments,
  and returns them as columns named \code{x0}, \code{y0}, \code{x1}
  and \code{y1} in a data frame. If the line segments were marked,
  the marks are appended as an extra column or columns to the
  data frame which is returned.  If the marks are a vector then a
  single column named \code{marks} is appended. in the data frame,
  with the same type as in the line segment pattern dataset.  If the
  marks are a data frame, then the columns of this data frame are
  appended (retaining their names).

}
\value{
  A data frame with 4 or 5 columns.
}
\examples{
  data(copper)
  df <- as.data.frame(copper$Lines)
}
\author{Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  and Rolf Turner \email{r.turner@auckland.ac.nz}
  
}
\keyword{spatial}
\keyword{manip}
