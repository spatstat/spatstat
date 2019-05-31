\name{psp}
\alias{psp}
\title{Create a Line Segment Pattern}
\description{
  Creates an object of class \code{"psp"} representing 
  a line segment pattern in the two-dimensional plane.
}
\usage{
  psp(x0,y0, x1, y1, window, marks=NULL,
        check=spatstat.options("checksegments"))
}
\arguments{
  \item{x0}{Vector of \eqn{x} coordinates of first endpoint of each segment}
  \item{y0}{Vector of \eqn{y} coordinates of first endpoint of each segment}
  \item{x1}{Vector of \eqn{x} coordinates of second endpoint of each segment}
  \item{y1}{Vector of \eqn{y} coordinates of second endpoint of each segment}
  \item{window}{window of observation,
    an object of class \code{"owin"}}
  \item{marks}{(optional) vector or data frame of mark values}
  \item{check}{Logical value indicating whether to check that the line segments
    lie inside the window.}
}
\value{
  An object of class \code{"psp"} 
  describing a line segment pattern in the two-dimensional plane
  (see \code{\link{psp.object}}).
}
\details{
  In the \pkg{spatstat} library, a spatial pattern of line segments is
  described by an object of class \code{"psp"}. This function
  creates such objects.

  The vectors \code{x0}, \code{y0}, \code{x1} and \code{y1} must be
  numeric vectors of equal length. They are interpreted as the cartesian
  coordinates of the endpoints of the line segments.

  A line segment pattern is assumed to have been observed within a specific
  region of the plane called the observation window.
  An object of class \code{"psp"} representing a point pattern
  contains information specifying the observation window.
  This window must always be specified when creating a point pattern dataset;
  there is intentionally no default action of ``guessing'' the window
  dimensions from the data points alone. 

  The argument \code{window} must be an object of class
  \code{"owin"}. It is a full description of the window geometry,
  and could have been obtained from \code{\link{owin}} or
  \code{\link{as.owin}}, or by just extracting the observation window
  of another dataset, or by manipulating such windows.
  See \code{\link{owin}} or the Examples below.

  The optional argument \code{marks} is given if the line segment pattern
  is marked, i.e. if each line segment carries additional information.
  For example, line segments which are classified into two or more different
  types, or colours, may be regarded as having a mark which identifies
  which colour they are. 

  The object \code{marks} must be a vector of the same length
  as \code{x0}, or a data frame with number of rows equal to the
  length of \code{x0}.  The interpretation is that \code{marks[i]}
  or \code{marks[i,]} is the mark attached to the \eqn{i}th line
  segment.  If the marks are real numbers then \code{marks} should
  be a numeric vector, while if the marks takes only a finite number
  of possible values (e.g. colours or types) then \code{marks}
  should be a \code{factor}.

  See \code{\link{psp.object}} for a description of the class
  \code{"psp"}.

  Users would normally invoke \code{psp} to create a line segment pattern,
  and the function \code{\link{as.psp}} to convert data in another
  format into a line segment pattern.
}
\seealso{
  \code{\link{psp.object}},
  \code{\link{as.psp}},
  \code{\link{owin.object}},
  \code{\link{owin}},
  \code{\link{as.owin}}.

  Function for extracting information from a segment pattern:
  \code{\link{marks.psp}},
  \code{\link{summary.psp}},
  \code{\link{midpoints.psp}},
  \code{\link{lengths.psp}}
  \code{\link{angles.psp}},
  \code{\link{endpoints.psp}}

  Convert line segments to infinite lines: \code{\link{extrapolate.psp}}.
}
\examples{
  X <- psp(runif(10), runif(10), runif(10), runif(10), window=owin())
  m <- data.frame(A=1:10, B=letters[1:10])
  X <- psp(runif(10), runif(10), runif(10), runif(10), window=owin(), marks=m)
}
\author{\adrian
  and \rolf.
}
\keyword{spatial}
\keyword{datagen}
