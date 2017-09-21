\name{clickjoin}
\alias{clickjoin}
\title{
  Interactively join vertices on a plot
}
\description{
  Given a point pattern representing a set of vertices,
  this command gives a point-and-click interface
  allowing the user to join pairs of selected vertices by edges.
}
\usage{
  clickjoin(X, \dots, add = TRUE, m = NULL, join = TRUE)
}
\arguments{
  \item{X}{
    Point pattern of vertices. An object of class \code{"ppp"}.
  }
  \item{\dots}{
    Arguments passed to \code{\link{segments}} to control
    the plotting of the new edges.
  }
  \item{add}{
    Logical. Whether the point pattern \code{X} should be added to the
    existing plot (\code{add=TRUE}) or a new plot should be created
    (\code{add=FALSE}).
  }
  \item{m}{
    Optional. Logical matrix specifying an initial
    set of edges. There is an edge between vertices \code{i} and
    \code{j} if \code{m[i,j] = TRUE}.
  }
  \item{join}{
    Optional. If \code{TRUE}, then each user click will join a pair of
    vertices. If \code{FALSE}, then each user click will delete an
    existing edge. This is only relevant if \code{m} is supplied.
  }
}
\details{
  This function makes it easier for the user to create a
  linear network or a planar graph, given a set of vertices.
  
  The function first displays the point pattern \code{X},
  then repeatedly prompts the user to click on a pair of points in \code{X}.
  Each selected pair of points will be joined
  by an edge. The function returns a
  logical matrix which has entries equal to \code{TRUE} for each
  pair of vertices joined by an edge.

  The selection of points is performed
  using \code{\link{identify.ppp}} which typically expects the user to
  click the left mouse button.  This point-and-click interaction
  continues until the user terminates it, 
  by pressing the middle mouse button,
  or pressing the right mouse button and selecting \code{stop}.

  The return value can be used in \code{\link{linnet}}
  to create a linear network.
}
\value{
  Logical matrix \code{m} with value \code{m[i,j] = TRUE} for every
  pair of vertices \code{X[i]} and \code{X[j]} that should be joined by
  an edge.
}
\author{
  \adrian.
}
\seealso{
  \code{\link{linnet}},
  \code{\link{clickppp}}
}
\keyword{spatial}
\keyword{datagen}
