\name{psib}
\alias{psib}
\alias{psib.kppm}
\title{
  Sibling Probability of Cluster Point Process
}
\description{
  Computes the sibling probability of a cluster point process model.
}
\usage{
  psib(object)

  \method{psib}{kppm}(object)
}
\arguments{
  \item{object}{
    Fitted cluster point process model
    (object of class \code{"kppm"}).
  }
}
\details{
  In a Poisson cluster process, two points are called \emph{siblings}
  if they belong to the same cluster, that is, if they had the same
  parent point. If two points of the process are
  separated by a distance \eqn{r}, the probability that
  they are siblings is \eqn{p(r) = 1 - 1/g(r)} where \eqn{g} is the
  pair correlation function of the process.
  
  The value \eqn{p(0) = 1 - 1/g(0)} is the probability that,
  if two points of the process are situated very close to each other,
  they came from the same cluster. This probability
  is an index of the strength of clustering, with high values
  suggesting strong clustering.

  This concept was proposed in Baddeley, Rubak and Turner (2015, page 479)
  and Baddeley (2016).
}
\value{
  A single number.
}
\references{
  Baddeley, A. (2016)
  Local composite likelihood for spatial point processes.
  \emph{Spatial Statistics}, in press.

  Baddeley, A., Rubak, E. and Turner, R. (2015)
  \emph{Spatial Point Patterns: Methodology and Applications with R}.
  Chapman and Hall/CRC Press.
}
\author{
  \adrian.
}
\seealso{
  \code{\link[spatstat]{kppm}}
}
\examples{
  fit <- kppm(redwood ~1, "Thomas")
  psib(fit)
}
\keyword{spatial}
\keyword{models}
