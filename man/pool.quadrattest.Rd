\name{pool.quadrattest}
\alias{pool.quadrattest}
\title{
  Pool Several Quadrat Tests
}
\description{
  Pool several quadrat tests into a single quadrat test.
}
\usage{
\method{pool}{quadrattest}(..., df=NULL, df.est=NULL, nsim=1999,
                                Xname=NULL, CR=NULL)
}
\arguments{
  \item{\dots}{
    Any number of objects, each of which is a quadrat test
    (object of class \code{"quadrattest"}).
  }
  \item{df}{
    Optional. Number of degrees of freedom of the test statistic.
    Relevant only for \eqn{\chi^2}{chi^2} tests.
    Incompatible with \code{df.est}.
  }
  \item{df.est}{
    Optional. The number of fitted parameters, or the
    number of degrees of freedom lost by estimation of
    parameters. 
    Relevant only for \eqn{\chi^2}{chi^2} tests.
    Incompatible with \code{df}.
  }
  \item{nsim}{
    Number of simulations, for Monte Carlo test.
  }
  \item{Xname}{
    Optional. Name of the original data.
  }
  \item{CR}{
    Optional. Numeric value of the Cressie-Read exponent \code{CR}
    overriding the value used in the tests.
  }
}
\details{
  The function \code{\link{pool}} is generic. This is the method for the
  class \code{"quadrattest"}.

  An object of class \code{"quadrattest"} represents a
  \eqn{\chi^2}{chi^2} test or Monte Carlo test
  of goodness-of-fit for a point process model, based on quadrat counts.
  Such objects are created by the command \code{\link{quadrat.test}}.

  Each of the arguments \code{\dots} must be an object of class
  \code{"quadrattest"}. They must all be the same type of test
  (chi-squared test or Monte Carlo test, conditional or unconditional)
  and must all have the same type of alternative hypothesis.

  The test statistic of the pooled test is the Pearson \eqn{X^2}
  statistic taken over all cells (quadrats) of all tests.
  The \eqn{p} value of the pooled test is then computed using
  either a Monte Carlo test or a \eqn{\chi^2}{chi^2} test.

  For a pooled \eqn{\chi^2}{chi^2} test, the number of degrees of freedom of
  the combined test is computed by adding the degrees of freedom
  of all the tests (equivalent to assuming the tests are independent)
  unless it is determined by the arguments \code{df} or \code{df.est}.
  The resulting \eqn{p} value is computed to obtain the
  pooled test.

  For a pooled Monte Carlo test, new simulations are performed
  to determine the pooled Monte Carlo \eqn{p} value. 
}
\value{
  Another object of class \code{"quadrattest"}.
}
\seealso{
  \code{\link{pool}},
  \code{\link{quadrat.test}}
}
\examples{
  Y <- split(humberside)
  test1 <- quadrat.test(Y[[1]])
  test2 <- quadrat.test(Y[[2]])
  pool(test1, test2, Xname="Humberside")
}
\author{Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  and Rolf Turner \email{r.turner@auckland.ac.nz}
  
}
\keyword{spatial}
\keyword{htest}

