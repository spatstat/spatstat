\name{heatkernelapprox}
\alias{heatkernelapprox}
\title{
  Approximation to Heat Kernel on Linear Network at Source Point
}
\description{
  Computes an approximation to the value of the heat kernel on a network
  evaluated at its source location.
}
\usage{
heatkernelapprox(X, sigma, nmax = 20, floored=TRUE)
}
\arguments{
  \item{X}{
    Point pattern on a linear network (object of class \code{"lpp"}).
  }
  \item{sigma}{
    Numeric. Bandwidth for kernel.
  }
  \item{nmax}{
    Number of terms to be used in the sum.
  }
  \item{floored}{
    Logical. If \code{TRUE}, all values are constrained to be
    greater than or equal to \eqn{1/L} where \eqn{L} is the total
    length of the network. This the exact value of the heat kernel
    when the bandwidth is infinite.
  }
}
\details{
  For each point \code{X[i]} in the pattern \code{X}, this algorithm computes an
  approximation to the value of the heat kernel with source point
  \code{X[i]} evaluated at the same location.
  
  The heat kernel \eqn{\kappa(u,v)} for a source location \eqn{u}
  evaluated at location \eqn{v}
  can be expressed as an infinite sum of contributions
  from all possible paths from \eqn{u} to \eqn{v}.
  This algorithm applies to the special case \eqn{u=v} where the source
  point and the query point are the same.

  The algorithm computes an approximation to \eqn{\kappa(u,u)}
  by taking only the contributions from paths 
  which (a) remain in the line segment containing the point \eqn{u} and
  (b) visit a vertex at most \code{nmax} times. 
}
\value{
  Numeric vector with one entry for each point in \code{X}.
}
\author{
  Greg McSwiggan and \adrian.
}
\seealso{
  \code{\link{hotrod}}
}
\examples{
  X <- runiflpp(3,simplenet)
  heatkernelapprox(X, 0.5)
}
\keyword{math}
