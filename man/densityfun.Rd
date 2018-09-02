\name{densityfun.ppp}
\alias{densityfun}
\alias{densityfun.ppp}
\title{
  Kernel Estimate of Intensity as a Spatial Function
}
\description{
  Compute a kernel estimate of intensity for a point pattern,
  and return the result as a function of spatial location.  
}
\usage{
densityfun(X, \dots)

\method{densityfun}{ppp}(X, sigma = NULL, \dots,
                        weights = NULL, edge = TRUE, diggle = FALSE)
}
\arguments{
  \item{X}{
    Point pattern (object of class \code{"ppp"}).
  }
  \item{sigma}{
    Smoothing bandwidth, or bandwidth selection function,
    passed to \code{\link{density.ppp}}.
  }
  \item{\dots}{
    Additional arguments passed to \code{\link{density.ppp}}.
  }
  \item{weights}{
    Optional vector of weights associated with the points of \code{X}.
  }
  \item{edge,diggle}{
    Logical arguments controlling the edge correction.
    Arguments passed to \code{\link{density.ppp}}.
  }
}
\details{
  The commands \code{densityfun} and \code{\link{density}}
  both perform kernel estimation of the intensity of a point pattern.
  The difference is that \code{\link{density}} returns a pixel image,
  containing the estimated intensity values at a grid of locations, while
  \code{densityfun} returns a \code{function(x,y)} which can be used
  to compute the intensity estimate at \emph{any} spatial location.
  For purposes such as model-fitting it is more accurate to
  use \code{densityfun}.
}
\value{
  A \code{function} with arguments \code{x,y}
  returning values of the intensity.
  The function also belongs to the class \code{"densityfun"} which has
  methods for \code{print} and \code{\link{as.im}}.
  It also belongs to the class \code{"funxy"} which has methods
  for \code{plot}, \code{contour} and \code{persp}.
}
\seealso{
  \code{\link{density}}.

  To interpolate values observed at the points, use \code{\link{Smoothfun}}.
}
\examples{
  f <- densityfun(swedishpines)
  f
  f(42, 60)
  plot(f)
}
\author{
  \spatstatAuthors.
}
\keyword{spatial}
\keyword{methods}
\keyword{smooth}
