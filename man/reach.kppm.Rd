\name{reach.kppm}
\alias{reach.kppm}
\title{Range of Interaction for a Cox or Cluster Point Process Model}
\description{
  Returns the range of interaction for a Cox or cluster point process model.
}
\usage{
  \method{reach}{kppm}(x, \dots, epsilon)
}
\arguments{
  \item{x}{
    Fitted point process model of class \code{"kppm"}.
  }
  \item{epsilon}{
    Optional numerical value. Differences smaller than
    \code{epsilon} are treated as zero.
  }
  \item{\dots}{
    Additional arguments passed to the range function
    of the given model.
  }
}
\details{
  The range of interaction for a fitted point process model
  of class \code{"kppm"}
  may defined as the smallest number \eqn{R} such that
  \eqn{g(r)=1} for all \eqn{r\ge R}{r>=R}, where \eqn{g} is the pair
  correlation function. 

  For many models the range is infinite, but one
  may instead use a value where the pair correlation function is
  sufficiently close to 1. The argument \code{epsilon} specifies the
  tolerance; there is a sensible default.
}
\value{Numeric}
\author{
  \spatstatAuthors.
}
\examples{
  fit <- kppm(redwood ~ 1)
  reach(fit)
}
\keyword{spatial}
\keyword{models}
