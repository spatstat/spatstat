\name{kernel.moment}
\alias{kernel.moment}
\title{Moment of Smoothing Kernel}
\description{
  Computes the complete or incomplete \eqn{m}th moment of a
  smoothing kernel.
}
\usage{
  kernel.moment(m, r, kernel = "gaussian")
}
\arguments{
  \item{m}{
    Exponent (order of moment).
    An integer.
  }
  \item{r}{
    Upper limit of integration for the incomplete moment.
    A numeric value or numeric vector.
    Set \code{r=Inf} to obtain the complete moment.
  }
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
  For more information about these kernels,
  see \code{\link[stats]{density.default}}. 

  The function \code{kernel.moment} computes the partial integral 
  \deqn{
    \int_{-\infty}^r t^m k(t) dt
  }{
    integral[-Inf][r] t^m k(t) dt
  }
  where \eqn{k(t)} is the selected kernel, \eqn{r} is the upper limit of
  integration, and \eqn{m} is the exponent or order.
  Here \eqn{k(t)} is the \bold{standard form} of the kernel,
  which has support \eqn{[-1,1]} and
  standard deviation \eqn{sigma = 1/c} where \code{c = kernel.factor(kernel)}.
}
\value{
  A single number, or a numeric vector of the same length as \code{r}.
}
\seealso{
  \code{\link[stats]{density.default}},
  \code{\link{dkernel}},
  \code{\link{kernel.factor}},
}
\examples{
   kernel.moment(1, 0.1, "epa")
   curve(kernel.moment(2, x, "epa"), from=-1, to=1)
}
\author{
  \adrian
  and Martin Hazelton.
}
\keyword{methods}
\keyword{nonparametric}
\keyword{smooth}
