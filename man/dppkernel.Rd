\name{dppkernel}
\alias{dppkernel}
\title{Extract Kernel from Determinantal Point Process Model Object}
\description{
  Returns the kernel of a determinantal point process model as a
  function of one argument \code{x}.
}
\usage{dppkernel(model, \dots)}
\arguments{
  \item{model}{Model of class \code{"detpointprocfamily"}.}
  \item{\dots}{Arguments passed to \code{\link{dppapproxkernel}} if the exact kernel is unknown}
}
\value{A function}
\author{
  \adrian
  
  
  \rolf
  
  and \ege
  
}
\examples{
kernelMatern <- dppkernel(dppMatern(lambda = 100, alpha=.01, nu=1, d=2))
plot(kernelMatern, xlim = c(0,0.1))
}
\keyword{spatial}
\keyword{models}
