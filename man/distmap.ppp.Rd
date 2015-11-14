\name{distmap.ppp}  %DontDeclareMethods
\alias{distmap.ppp}
\title{
  Distance Map of Point Pattern
}
\description{
  Computes the distance from each pixel to the nearest
  point in the given point pattern.
}
\usage{
  \method{distmap}{ppp}(X, \dots)
}
\arguments{
  \item{X}{A point pattern (object of class \code{"ppp"}).
  }
  \item{\dots}{Arguments passed to \code{\link{as.mask}}
    to control pixel resolution.
  }
}
\value{
  A pixel image (object of class \code{"im"}) whose greyscale values
  are the values of the distance map.
  The return value has attributes \code{"index"} and \code{"bdry"}
  which are also pixel images. 
}
\details{
  The ``distance map'' of a point pattern \eqn{X} is the function
  \eqn{f} whose value \code{f(u)} is defined for any two-dimensional
  location \eqn{u} as the shortest distance from \eqn{u} to \eqn{X}.

  This function computes the distance map of the point pattern \code{X}
  and returns the distance map as a pixel image. The greyscale value
  at a pixel \eqn{u} equals the distance from \eqn{u} to the nearest
  point of the pattern \code{X}.

  Additionally, the return value 
  has two attributes, \code{"index"} and \code{"bdry"}, which are
  also pixel images. The grey values in \code{"bdry"} give the
  distance from each pixel to the bounding rectangle of the image.
  The grey values in \code{"index"} are integers identifying which
  point of \code{X} is closest. 
  
  This is a method for the generic function \code{\link{distmap}}.
  
  Note that this function gives the distance from the
  \emph{centre of each pixel} to the nearest data point.
  To compute the exact distance from a given spatial location
  to the nearest data point in \code{X}, use \code{\link{distfun}} or 
  \code{\link{nncross}}.
}
\seealso{
  Generic function \code{\link{distmap}} and other methods
  \code{\link{distmap.psp}},
  \code{\link{distmap.owin}}.
  
  Generic function \code{\link{distfun}}.

  Nearest neighbour distance \code{\link{nncross}}
}
\examples{
  data(cells)
  U <- distmap(cells)
  \dontrun{
  plot(U)
  plot(attr(U, "bdry"))
  plot(attr(U, "index"))
  }
}
\author{\adrian
  
  
  and \rolf
  
}
\keyword{spatial}
\keyword{math}
