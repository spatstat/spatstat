\name{kernel.factor}
\alias{kernel.factor}
\title{Scale factor for density kernel} 
\description{
  Returns a scale factor for the kernels used in density estimation
  for numerical data.
}
\usage{
  kernel.factor(kernel = "gaussian")
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
}
\details{
  Kernel estimation of a probability density in one dimension
  is performed by \code{\link[stats]{density.default}}
  using a kernel function selected from the list above.

  This function computes a scale constant for the kernel.
  For the Gaussian kernel, this constant is equal to 1.
  Otherwise, the constant \eqn{c} is such that the kernel
  with standard deviation \eqn{1} is supported on the interval
  \eqn{[-c,c]}. 
  
  For more information about these kernels,
  see \code{\link[stats]{density.default}}. 
}
\value{
  A single number.
}
\seealso{
  \code{\link[stats]{density.default}},
  \code{\link{dkernel}},
  \code{\link{kernel.moment}},
  \code{\link{kernel.squint}}
}
\examples{
   kernel.factor("rect")
   # bandwidth for Epanechnikov kernel with half-width h=1
   h <- 1
   bw <- h/kernel.factor("epa")
}
\author{\adrian
  and Martin Hazelton
}
\keyword{methods}
\keyword{nonparametric}
\keyword{smooth}
