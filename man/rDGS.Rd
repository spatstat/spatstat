\name{rDGS}
\alias{rDGS}
\title{Perfect Simulation of the Diggle-Gates-Stibbard Process}
\description{
  Generate a random pattern of points, a simulated realisation
  of the Diggle-Gates-Stibbard process, using a perfect simulation algorithm.
}
\usage{
  rDGS(beta, rho, W = owin(), expand=TRUE, nsim=1, drop=TRUE)
}
\arguments{
  \item{beta}{
    intensity parameter (a positive number).
  }
  \item{rho}{
    interaction range (a non-negative number).
  }
  \item{W}{
    window (object of class \code{"owin"}) in which to
    generate the random pattern. 
  }
  \item{expand}{
    Logical. If \code{FALSE}, simulation is performed
    in the window \code{W}, which must be rectangular.
    If \code{TRUE} (the default), simulation is performed
    on a larger window, and the result is clipped to the original
    window \code{W}.
    Alternatively \code{expand} can be an object of class 
    \code{"rmhexpand"} (see \code{\link{rmhexpand}})
    determining the expansion method.
  }
  \item{nsim}{Number of simulated realisations to be generated.}  
  \item{drop}{
    Logical. If \code{nsim=1} and \code{drop=TRUE} (the default), the
    result will be a point pattern, rather than a list 
    containing a point pattern.
  }
}
\details{
  This function generates a realisation of the 
  Diggle-Gates-Stibbard point process in the window \code{W}
  using a \sQuote{perfect simulation} algorithm.

  Diggle, Gates and Stibbard (1987) proposed a pairwise interaction
  point process in which each pair of points separated by
  a distance \eqn{d} contributes a factor \eqn{e(d)} to the
  probability density, where
  \deqn{
    e(d) = \sin^2\left(\frac{\pi d}{2\rho}\right)
  }{
    e(d) = sin^2((pi * d)/(2 * rho))
  }
  for \eqn{d < \rho}{d < rho}, and \eqn{e(d)} is equal to 1
  for \eqn{d \ge \rho}{d >= rho}.

  The simulation algorithm used to generate the point pattern
  is \sQuote{dominated coupling from the past}
  as implemented by Berthelsen and \Moller (2002, 2003).
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
  If \code{nsim = 1}, a point pattern (object of class \code{"ppp"}).
  If \code{nsim > 1}, a list of point patterns.
}
\references{
  Berthelsen, K.K. and \Moller, J. (2002)
  A primer on perfect simulation for spatial point processes.
  \emph{Bulletin of the Brazilian Mathematical Society} 33, 351-367.

  Berthelsen, K.K. and \Moller, J. (2003)
  Likelihood and non-parametric Bayesian MCMC inference 
  for spatial point processes based on perfect simulation and
  path sampling. 
  \emph{Scandinavian Journal of Statistics} 30, 549-564.

  Diggle, P.J., Gates, D.J., and Stibbard, A. (1987)
   A nonparametric estimator for pairwise-interaction point processes.
   Biometrika \bold{74}, 763 -- 770.
   \emph{Scandinavian Journal of Statistics} \bold{21}, 359--373.

  \Moller, J. and Waagepetersen, R. (2003).
  \emph{Statistical Inference and Simulation for Spatial Point Processes.}
  Chapman and Hall/CRC.
}
\author{
  \adrian,
  based on original code for the Strauss process by 
  Kasper Klitgaard Berthelsen.
}
\examples{
   X <- rDGS(50, 0.05)
}
\seealso{
  \code{\link{rmh}},
  \code{\link{DiggleGatesStibbard}}.

  \code{\link{rStrauss}},
  \code{\link{rHardcore}},
  \code{\link{rStraussHard}},
  \code{\link{rDiggleGratton}},
  \code{\link{rPenttinen}}.
}
\keyword{spatial}
\keyword{datagen}
