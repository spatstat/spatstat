\name{rDiggleGratton}
\alias{rDiggleGratton}
\title{Perfect Simulation of the Diggle-Gratton Process}
\description{
  Generate a random pattern of points, a simulated realisation
  of the Diggle-Gratton process, using a perfect simulation algorithm.
}
\usage{
  rDiggleGratton(beta, delta, rho, kappa=1, W = owin())
}
\arguments{
  \item{beta}{
    intensity parameter (a positive number).
  }
  \item{delta}{
    hard core distance (a non-negative number).
  }
  \item{rho}{
    interaction range (a number greater than \code{delta}).
  }
  \item{kappa}{
    interaction exponent (a non-negative number).
  }
  \item{W}{
    window (object of class \code{"owin"}) in which to
    generate the random pattern. Currently this must be a rectangular
    window.
  }
}
\details{
  This function generates a realisation of the 
  Diggle-Gratton point process in the window \code{W}
  using a \sQuote{perfect simulation} algorithm.

  Diggle and Gratton (1984, pages 208-210)
  introduced the pairwise interaction point
  process with pair potential \eqn{h(t)} of the form
  \deqn{
    h(t) = \left( \frac{t-\delta}{\rho-\delta} \right)^\kappa
    \quad\quad \mbox{  if  } \delta \le t \le \rho
  }{
    h(t) = ((t - delta)/(rho - delta))^kappa, {    } delta <= t <= rho
  }
  with \eqn{h(t) = 0} for \eqn{t < \delta}{t < delta}
  and  \eqn{h(t) = 1} for \eqn{t > \rho}{t > rho}.
  Here \eqn{\delta}{delta}, \eqn{\rho}{rho} and \eqn{\kappa}{kappa}
  are parameters.

  Note that we use the symbol \eqn{\kappa}{kappa}
  where Diggle and Gratton (1984) 
  use \eqn{\beta}{beta}, since in \pkg{spatstat} we reserve the symbol
  \eqn{\beta}{beta} for an intensity parameter.

  The parameters must all be nonnegative,
  and must satisfy \eqn{\delta \le \rho}{delta <= rho}.

  The simulation algorithm used to generate the point pattern
  is \sQuote{dominated coupling from the past}
  as implemented by Berthelsen and Moller (2002, 2003).
  This is a \sQuote{perfect simulation} or \sQuote{exact simulation}
  algorithm, so called because the output of the algorithm is guaranteed
  to have the correct probability distribution exactly (unlike the
  Metropolis-Hastings algorithm used in \code{\link{rmh}}, whose output
  is only approximately correct).

  There is a tiny chance that the algorithm will 
  run out of space before it has terminated. If this occurs, an error
  message will be generated.
}
\value{
  A point pattern (object of class \code{"ppp"}).
}
\references{
  Berthelsen, K.K. and Moller, J. (2002)
  A primer on perfect simulation for spatial point processes.
  \emph{Bulletin of the Brazilian Mathematical Society} 33, 351-367.

  Berthelsen, K.K. and Moller, J. (2003)
  Likelihood and non-parametric Bayesian MCMC inference 
  for spatial point processes based on perfect simulation and
  path sampling. 
  \emph{Scandinavian Journal of Statistics} 30, 549-564.

  Diggle, P.J. and Gratton, R.J. (1984)
  Monte Carlo methods of inference for implicit statistical models.
  \emph{Journal of the Royal Statistical Society, series B}
  \bold{46}, 193 -- 212.

  Moller, J. and Waagepetersen, R. (2003).
  \emph{Statistical Inference and Simulation for Spatial Point Processes.}
  Chapman and Hall/CRC.
}
\author{
  Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  based on original code for the Strauss process by 
  Kasper Klitgaard Berthelsen.
}
\examples{
   X <- rDiggleGratton(50, 0.02, 0.07)
}
\seealso{
  \code{\link{rmh}},
  \code{\link{DiggleGratton}},
  \code{\link{rStrauss}},
  \code{\link{rHardcore}}.
}
\keyword{spatial}
\keyword{datagen}
