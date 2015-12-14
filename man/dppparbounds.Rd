\name{dppparbounds}
\alias{dppparbounds}
\title{Parameter Bound for a Determinantal Point Process Model}
\description{
  Returns the lower and upper bound for a specific parameter of a
  determinantal point process model when all other parameters are
  fixed.
}
\usage{dppparbounds(model, name, \dots)}
\arguments{
  \item{model}{Model of class \code{"detpointprocfamily"}.}
  \item{name}{name of the parameter for which the bound should be computed.}
  \item{\dots}{
    Additional arguments passed to the \code{parbounds} function of the given
    model
  }
}
\value{A \code{data.frame} containing lower and upper bounds.}
\author{
  Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  Rolf Turner \email{r.turner@auckland.ac.nz}
  
  and Ege Rubak \email{rubak@math.aau.dk}
  
}
\examples{
model <- dppMatern(lambda=100, alpha=.01, nu=1, d=2)
dppparbounds(model, "lambda")
}
\keyword{spatial}
\keyword{models}
