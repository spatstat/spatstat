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
  \code{\link{dkernel}}
}
\examples{
   kernel.factor("epa")
}
\author{Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  \email{adrian@maths.uwa.edu.au}
  
  and Martin Hazelton
}
\keyword{methods}
\keyword{nonparametric}
\keyword{smooth}
