\name{psp.object}
\alias{psp.object} %DoNotExport
\title{Class of Line Segment Patterns}
\description{
  A class \code{"psp"} to represent a spatial pattern of
  line segments in the plane.
  Includes information about the window in which the
  pattern was observed. Optionally includes marks.
}
\details{
  An object of this class represents 
  a two-dimensional pattern of line segments. It specifies
  \itemize{
    \item the locations of the line segments (both endpoints)
    \item the window in which the pattern was observed
    \item optionally, a ``mark'' attached to each line segment
      (extra information such as a type label).
    }
  If \code{X} is an object of type \code{psp},
  it contains the following elements:
  \tabular{ll}{
    \code{ends} \tab data frame with entries \code{x0, y0, x1, y1} \cr
                \tab giving coordinates of segment endpoints \cr 
    \code{window} \tab window of observation \cr
                  \tab (an object of class \code{\link{owin}}) \cr
    \code{n} \tab number of line segments \cr
    \code{marks} \tab optional vector or data frame of marks \cr
    \code{markformat} \tab character string specifying the format of the \cr
                      \tab marks; \dQuote{none}, \dQuote{vector}, or
                      \dQuote{dataframe}
  }
  Users are strongly advised not to manipulate these entries
  directly.
  
  Objects of class \code{"psp"}
  may be created by the function
  \code{\link{psp}}
  and converted from other types of data by the function
  \code{\link{as.psp}}.
  Note that you must always specify the window of observation;
  there is intentionally no default action of ``guessing'' the window
  dimensions from the line segments alone.

  Subsets of a line segment pattern may be obtained by the functions
  \code{\link{[.psp}} and \code{\link{clip.psp}}.

  Line segment pattern objects can be plotted just by typing \code{plot(X)}
  which invokes the \code{plot} method for line segment pattern objects,
  \code{\link{plot.psp}}. See \code{\link{plot.psp}} for further information.

  There are also methods for \code{summary} and \code{print}
  for line segment patterns. Use \code{summary(X)} to see a useful description
  of the data.

  Utilities for line segment patterns include
  \code{\link{midpoints.psp}} (to compute the midpoints of each segment),
  \code{\link{lengths_psp}}, (to compute the length of each segment),
  \code{\link{angles.psp}}, (to compute the angle of orientation of
  each segment), and 
  \code{\link{distmap.psp}} to compute the distance map of a
  line segment pattern.
}
\seealso{
  \code{\link{psp}},
  \code{\link{as.psp}},
  \code{\link{[.psp}}
}
\examples{
# creating 
    a <- psp(runif(20),runif(20),runif(20),runif(20), window=owin())
# converting from other formats
    a <- as.psp(matrix(runif(80), ncol=4), window=owin())
    a <- as.psp(data.frame(x0=runif(20), y0=runif(20),
                            x1=runif(20), y1=runif(20)), window=owin())
# clipping
    w <- owin(c(0.1,0.7), c(0.2, 0.8))
    b <- clip.psp(a, w)
    b <- a[w]
# the last two lines are equivalent.
}
\author{\adrian
  
  
  and \rolf
  
}
\keyword{spatial}
\keyword{attribute}
 
 
