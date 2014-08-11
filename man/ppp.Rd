\name{ppp}
\alias{ppp}
\title{Create a Point Pattern}
\description{
  Creates an object of class \code{"ppp"} representing 
  a point pattern dataset in the two-dimensional plane.
}
\usage{
  ppp(x,y, \dots, window, marks, check=TRUE)
}
\arguments{
  \item{x}{Vector of \eqn{x} coordinates of data points}
  \item{y}{Vector of \eqn{y} coordinates of data points}
  \item{window}{window of observation,
    an object of class \code{"owin"}}
  \item{\dots}{arguments passed to \code{\link{owin}} to create the
    window, if \code{window} is missing}
  \item{marks}{(optional) mark values for the points.
    A vector or data frame.}
  \item{check}{Logical flag indicating whether to check
    that all the \eqn{(x,y)} points lie inside the specified window.
    Do not set this to \code{FALSE} unless you are sure that this
    check is unnecessary.
  }
}
\value{
  An object of class \code{"ppp"} 
  describing a point pattern in the two-dimensional plane
  (see \code{\link{ppp.object}}).
}
\details{
  In the \pkg{spatstat} library, a point pattern dataset is
  described by an object of class \code{"ppp"}. This function
  creates such objects.

  The vectors \code{x} and \code{y} must be numeric vectors of
  equal length. They are interpreted as the cartesian coordinates
  of the points in the pattern.

  A point pattern dataset is assumed to have been observed within a specific
  region of the plane called the observation window.
  An object of class \code{"ppp"} representing a point pattern
  contains information specifying the observation window.
  This window must always be specified when creating a point pattern dataset;
  there is intentionally no default action of ``guessing'' the window
  dimensions from the data points alone. 

  You can specify the observation window in several
  (mutually exclusive) ways:
  \itemize{
    \item
    \code{xrange, yrange} specify a rectangle
    with these dimensions;
    \item
    \code{poly} specifies a polygonal boundary.
    If the boundary is a single polygon then \code{poly}
    must be a list with components \code{x,y}
    giving the coordinates of the vertices.
    If the boundary consists of several disjoint polygons
    then \code{poly} must be a list of such lists
    so that \code{poly[[i]]$x} gives the \eqn{x} coordinates
    of the vertices of the \eqn{i}th boundary polygon.
    \item
    \code{mask} specifies a binary pixel image with entries
    that are \code{TRUE} if the corresponding pixel is inside
    the window.
    \item
    \code{window} is an object of class \code{"owin"}
    (see \code{\link{owin.object}}) specifying the window.
  }
  The arguments \code{xrange, yrange} or \code{poly}
  or \code{mask} are passed to the window creator function
  \code{\link{owin}} for interpretation. See
  \code{\link{owin}} for further details.

  The argument \code{window}, if given, must be an object of class
  \code{"owin"}. It is a full description of the window geometry,
  and could have been obtained from \code{\link{owin}} or
  \code{\link{as.owin}}, or by just extracting the observation window
  of another point pattern, or by manipulating such windows.
  See \code{\link{owin}} or the Examples below.

  The points with coordinates \code{x} and \code{y}
  \bold{must} lie inside the specified window, in order to
  define a valid object of this class. 
  Any points which do not lie inside the window will be
  removed from the point pattern, and a warning will be issued.
  See the section on Rejected Points.

  The name of the unit of length for the \code{x} and \code{y} coordinates
  can be specified in the dataset, using the argument \code{unitname}, which is
  passed to \code{\link{owin}}. See the examples below, or the help file
  for \code{\link{owin}}.
  
  The optional argument \code{marks} is given if the point pattern
  is marked, i.e. if each data point carries additional information.
  For example, points which are classified into two or more different
  types, or colours, may be regarded as having a mark which identifies
  which colour they are. Data recording the locations and heights of
  trees in a forest can be regarded as a marked point pattern where the
  mark is the tree height.
  
  The argument \code{marks} can be either
  \itemize{
    \item a vector, of
    the same length as \code{x} and \code{y}, which is interpreted so
    that \code{marks[i]} is the mark attached to the point
    \code{(x[i],y[i])}. If the mark is a real number then \code{marks}
    should be a numeric vector, while if the mark takes only a finite
    number of possible values (e.g. colours or types) then
    \code{marks} should be a \code{factor}.
    \item
    a data frame, with the number of rows equal to the number of points
    in the point pattern. The \code{i}th row of the data frame is interpreted
    as containing the mark values for the \code{i}th point in the point
    pattern. The columns of the data frame correspond to different
    mark variables (e.g. tree species and tree diameter).
  }
  
  See \code{\link{ppp.object}} for a description of the
  class \code{"ppp"}.

  Users would normally invoke \code{ppp} to create a point pattern,
  but the functions \code{\link{as.ppp}} and 
  \code{scanpp} may sometimes be convenient.
}
\section{Rejected points}{
  The points with coordinates \code{x} and \code{y}
  \bold{must} lie inside the specified window, in order to
  define a valid object of class \code{"ppp"}.
  Any points which do not lie inside the window will be
  removed from the point pattern, and a warning will be issued.

  The rejected points are still accessible: they are stored
  as an attribute of the point pattern called \code{"rejects"}
  (which is an object of class \code{"ppp"} containing the rejected points
  in a large window). However, rejected points in a point pattern
  will be ignored by all other functions except
  \code{\link{plot.ppp}}.

  To remove the rejected points altogether,
  use \code{\link{as.ppp}}. To include the rejected points,
  you will need to find a larger window that contains them,
  and use this larger window in a call to \code{ppp}.
}

\seealso{
  \code{\link{ppp.object}},
  \code{\link{as.ppp}},
  \code{\link{owin.object}},
  \code{\link{owin}},
  \code{\link{as.owin}}
}
\examples{
  # some arbitrary coordinates in [0,1]
  x <- runif(20)
  y <- runif(20)

  # the following are equivalent
  X <- ppp(x, y, c(0,1), c(0,1))
  X <- ppp(x, y)
  X <- ppp(x, y, window=owin(c(0,1),c(0,1)))

  # specify that the coordinates are given in metres
  X <- ppp(x, y, c(0,1), c(0,1), unitname=c("metre","metres"))

  \dontrun{plot(X)}

  # marks
  m <- sample(1:2, 20, replace=TRUE)
  m <- factor(m, levels=1:2)
  X <- ppp(x, y, c(0,1), c(0,1), marks=m)
  \dontrun{plot(X)}

  # polygonal window
  X <- ppp(x, y, poly=list(x=c(0,10,0), y=c(0,0,10)))
  \dontrun{plot(X)}

  # copy the window from another pattern
  data(cells)
  X <- ppp(x, y, window=cells$window)
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{datagen}
