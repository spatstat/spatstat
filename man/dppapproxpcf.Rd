\name{dppapproxpcf}
\alias{dppapproxpcf}
\title{Approximate Pair Correlation Function of Determinantal Point Process Model}
\description{
  Returns an approximation to the
  theoretical pair correlation function of a
  determinantal point process model, as a
  function of one argument \eqn{x}.
}
\usage{dppapproxpcf(model, trunc = 0.99, W = NULL)}
\arguments{
  \item{model}{Object of class \code{"detpointprocfamily"}.}
  \item{trunc}{Numeric specifying how the model truncation is
    performed. See Details section of \code{\link{simulate.detpointprocfamily}}.}
  \item{W}{Optional window -- undocumented at the moment.}
}
\details{This function is usually NOT needed for anything. It
only exists for investigative purposes.}
\author{
  Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  Rolf Turner \email{r.turner@auckland.ac.nz}
  
  and Ege Rubak \email{rubak@math.aau.dk}
  
}
\examples{
f <- dppapproxpcf(dppMatern(lambda = 100, alpha=.028, nu=1, d=2))
plot(f, xlim = c(0,0.1))
}
