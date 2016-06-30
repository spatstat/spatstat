\name{Smoothfun.ppp}
\alias{Smoothfun}
\alias{Smoothfun.ppp}
\title{
  Smooth Interpolation of Marks as a Spatial Function
}
\description{
  Perform spatial smoothing of numeric values observed
  at a set of irregular locations, and return the result
  as a function of spatial location.  
}
\usage{
Smoothfun(X, \dots)

\method{Smoothfun}{ppp}(X, sigma = NULL, \dots,
                        weights = NULL, edge = TRUE, diggle = FALSE)
}
\arguments{
  \item{X}{
    Marked point pattern (object of class \code{"ppp"}).
  }
  \item{sigma}{
    Smoothing bandwidth, or bandwidth selection function,
    passed to \code{\link{Smooth.ppp}}.
  }
  \item{\dots}{
    Additional arguments passed to \code{\link{Smooth.ppp}}.
  }
  \item{weights}{
    Optional vector of weights associated with the points of \code{X}.
  }
  \item{edge,diggle}{
    Logical arguments controlling the edge correction.
    Arguments passed to \code{\link{Smooth.ppp}}.
  }
}
\details{
  The commands \code{Smoothfun} and \code{\link{Smooth}}
  both perform kernel-smoothed spatial interpolation
  of numeric values observed at irregular spatial locations.
  The difference is that \code{\link{Smooth}} returns a pixel image,
  containing the interpolated values at a grid of locations, while
  \code{Smoothfun} returns a \code{function(x,y)} which can be used
  to compute the interpolated value at \emph{any} spatial location.
  For purposes such as model-fitting it is more accurate to
  use \code{Smoothfun} to interpolate data.
}
\value{
  A \code{function} with arguments \code{x,y}.
  The function also belongs to the class \code{"Smoothfun"} which has
  methods for \code{print} and \code{\link{as.im}}.
  It also belongs to the class \code{"funxy"} which has methods
  for \code{plot}, \code{contour} and \code{persp}.
}
\seealso{
  \code{\link{Smooth}}
}
\examples{
  f <- Smoothfun(longleaf)
  f
  f(120, 80)
  plot(f)
}
\author{\adrian
  
  ,
  \rolf
  
  and \ege.
}
\keyword{spatial}
\keyword{methods}
\keyword{smooth}
