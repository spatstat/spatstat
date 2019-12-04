\name{reach.dppm}
\alias{reach.dppm}
\alias{reach.detpointprocfamily}
\title{Range of Interaction for a Determinantal Point Process Model}
\description{
  Returns the range of interaction for a determinantal point process model.
}
\usage{
  \method{reach}{dppm}(x, \dots)

  \method{reach}{detpointprocfamily}(x, \dots)
}
\arguments{
  \item{x}{Model of class \code{"detpointprocfamily"} or \code{"dppm"}.}
  \item{\dots}{Additional arguments passed to the range function
    of the given model.
  }
}
\details{
  The range of interaction for a determinantal point
  process model may defined as the smallest number \eqn{R} such that
  \eqn{g(r)=1} for all \eqn{r\ge R}{r>=R}, where \eqn{g} is the pair
  correlation function. For many models the range is infinite, but one
  may instead use a value where the pair correlation function is
  sufficiently close to 1. For example in the \Matern model this defaults
  to finding \eqn{R} such that \eqn{g(R)=0.99}.
}
\value{Numeric}
\author{
  \spatstatAuthors.
}
\examples{
reach(dppMatern(lambda=100, alpha=.01, nu=1, d=2))
}
\keyword{spatial}
\keyword{models}
