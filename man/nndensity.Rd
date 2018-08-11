\name{nndensity.ppp}
\alias{nndensity}
\alias{nndensity.ppp}
\title{
  Estimate Intensity of Point Pattern Using Nearest Neighbour Distances
}
\description{
  Estimates the intensity of a point pattern
  using the distance from each spatial location
  to the \code{k}th nearest data point.
}
\usage{
nndensity(x, ...)

\method{nndensity}{ppp}(x, k, ..., verbose = TRUE)
}
\arguments{
  \item{x}{
    A point pattern (object of class \code{"ppp"})
    or some other spatial object.
  }
  \item{k}{
    Integer. The distance to the \code{k}th nearest data point
    will be computed. There is a sensible default.
  }
  \item{\dots}{
    Arguments passed to \code{\link{nnmap}} and \code{\link{as.mask}}
    controlling the pixel resolution.
  }
  \item{verbose}{
    Logical. If \code{TRUE}, print the value of \code{k}
    when it is automatically selected. If \code{FALSE}, remain silent.
  }
}
\details{
  This function computes a quick estimate of the intensity of the point
  process that generated the point pattern \code{x}.

  For each spatial location \eqn{s}, let \eqn{d(s)} be the distance from \eqn{s}
  to the \eqn{k}-th nearest point in the dataset \code{x}.
  If the data came from a homogeneous
  Poisson process with intensity \eqn{\lambda}{lambda},
  then \eqn{\pi d(s)^2}{pi * d(s)^2} would follow a 
  negative exponential distribution with mean
  \eqn{1/\lambda}{1/lambda}, and the maximum likelihood estimate of
  \eqn{\lambda}{lambda} would be \eqn{1/(\pi d(s)^2)}{1/(pi * d(s)^2)}.
  This is the estimate computed by \code{nndensity},
  apart from an edge effect correction.

  This estimator of intensity is relatively fast to compute, 
  and is spatially adaptive
  (so that it can handle wide variation in the intensity
  function). However, it implicitly assumes the points are independent,
  so it does not perform well if the pattern is strongly clustered
  or strongly inhibited.

  The value of \code{k} should be greater than 1 in order to avoid
  infinite peaks in the intensity estimate around each data point.
  The default value of \code{k} is the square root of the number of
  points in \code{x}, which seems to work well in many cases.
  
  The window of \code{x} is digitised using \code{\link{as.mask}}
  and the values \eqn{d(s)} are computed using \code{\link{nnmap}}.
  To control the pixel resolution, see \code{\link{as.mask}}.
}
\value{
  A pixel image (object of class \code{"im"}) giving the
  estimated intensity of the point process at each spatial location.
  Pixel values are intensities (number of points per unit area).
}
\references{
  NEED REFERENCES. TRY CRESSIE 
}
\seealso{
 \code{\link{density.ppp}},
 \code{\link{intensity}}
 for alternative estimates of point process intensity.
}
\examples{
   plot(nndensity(swedishpines))
}
\author{\adrian
  
  
  and \rolf
  
}
\keyword{spatial}
\keyword{methods}
\keyword{smooth}
