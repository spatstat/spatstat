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
  \item{trunc}{Numeric value specifying how the model truncation is
    performed.
    See Details section of \code{\link{simulate.detpointprocfamily}}.}
  \item{W}{Optional window -- undocumented at the moment.}
}
\details{
  This function is usually NOT needed for anything. It
  only exists for investigative purposes.
}
\value{
  A function in the \R language, with one numeric argument \code{x},
  that returns the value of the approximate pair correlation
  at distances \code{x}.
}
\author{
  \spatstatAuthors.
}
\examples{
f <- dppapproxpcf(dppMatern(lambda = 100, alpha=.028, nu=1, d=2))
plot(f, xlim = c(0,0.1))
}
