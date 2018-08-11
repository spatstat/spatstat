\name{kernel.squint}
\alias{kernel.squint}
\title{Integral of Squared Kernel} 
\description{
  Computes the integral of the squared kernel,
  for the kernels used in density estimation
  for numerical data.
}
\usage{
  kernel.squint(kernel = "gaussian", bw=1)
}
\arguments{
  \item{kernel}{
    String name of the kernel.
    Options are
    \code{"gaussian"}, \code{"rectangular"},
    \code{"triangular"},
    \code{"epanechnikov"},
    \code{"biweight"},
    \code{"cosine"} and \code{"optcosine"}.
    (Partial matching is used).
  }
  \item{bw}{
    Bandwidth (standard deviation) of the kernel.
  }
}
\details{
  Kernel estimation of a probability density in one dimension
  is performed by \code{\link[stats]{density.default}}
  using a kernel function selected from the list above.

  This function computes the integral of the squared kernel,
  \deqn{
    R = \int_{-\infty}^{\infty} k(x)^2 \, {\rm d}x
  }{
    R = integral of k(x)^2 dx from x = -infinity to x = +infinity
  }
  where \eqn{k(x)} is the kernel with bandwidth \code{bw}.
}
\value{
  A single number.
}
\seealso{
  \code{\link[stats]{density.default}},
  \code{\link{dkernel}},
  \code{\link{kernel.moment}},
  \code{\link{kernel.factor}}
}
\examples{
   kernel.squint("gaussian", 3)

   # integral of squared Epanechnikov kernel with half-width h=1
   h <- 1
   bw <- h/kernel.factor("epa")
   kernel.squint("epa", bw)
}
\author{
  \spatstatAuthors
  and Martin Hazelton
}
\keyword{methods}
\keyword{nonparametric}
\keyword{smooth}
