\name{circdensity}
\alias{circdensity}
\title{
  Density Estimation for Circular Data
}
\description{
  Computes a kernel smoothed estimate of the probability density
  for angular data.
}
\usage{
circdensity(x, sigma = "nrd0", \dots,
               bw = NULL,
               weights=NULL, unit = c("degree", "radian"))
}
\arguments{
  \item{x}{
    Numeric vector, containing angular data.
  }
  \item{sigma}{
    Smoothing bandwidth, or bandwidth selection rule, passed to
    \code{\link[stats]{density.default}}.
  }
  \item{bw}{Alternative to \code{sigma} for consistency with other functions.}
  \item{\dots}{
    Additional arguments passed to
    \code{\link[stats]{density.default}},
    such as \code{kernel} and \code{weights}.
  }
  \item{weights}{
    Optional numeric vector of weights for the data in \code{x}.
  }
  \item{unit}{
    The unit of angle in which \code{x} is expressed.
  }
}
\details{
  The angular values \code{x} are smoothed using
  (by default) the wrapped Gaussian kernel with standard deviation \code{sigma}.
}
\value{
  An object of class \code{"density"} (produced by
  \code{\link[stats]{density.default}}) which can be plotted
  by \code{plot} or by \code{\link{rose}}.
}
\author{\adrian
  
  
  \rolf
  
  and \ege
  
}
\seealso{
  \code{\link[stats]{density.default}}),
  \code{\link{rose}}.
}
\examples{
  ang <- runif(1000, max=360)
  rose(circdensity(ang, 12))
}
\keyword{nonparametric}
\keyword{smooth}
