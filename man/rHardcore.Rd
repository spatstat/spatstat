\name{rHardcore}
\alias{rHardcore}
\title{Perfect Simulation of the Hardcore Process}
\description{
  Generate a random pattern of points, a simulated realisation
  of the Hardcore process, using a perfect simulation algorithm.
}
\usage{
  rHardcore(beta, R = 0, W = owin(), expand=TRUE, nsim=1, drop=TRUE)
}
\arguments{
  \item{beta}{
    intensity parameter (a positive number).
  }
  \item{R}{
    hard core distance (a non-negative number).
  }
  \item{W}{
    window (object of class \code{"owin"}) in which to
    generate the random pattern. Currently this must be a rectangular
    window.
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
  Hardcore point process in the window \code{W}
  using a \sQuote{perfect simulation} algorithm.

  The Hardcore process is a model for strong spatial inhibition.
  Two points of the process are forbidden to lie closer than
  \code{R} units apart.
  The Hardcore process is the special case of the Strauss process
  (see \code{\link{rStrauss}})
  with interaction parameter \eqn{\gamma}{gamma} equal to zero.
  
  The simulation algorithm used to generate the point pattern
  is \sQuote{dominated coupling from the past}
  as implemented by Berthelsen and \ifelse{latex}{\out{M\o ller}}{Moller} (2002, 2003).
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
  Berthelsen, K.K. and \ifelse{latex}{\out{M\o ller}}{Moller}, J. (2002)
  A primer on perfect simulation for spatial point processes.
  \emph{Bulletin of the Brazilian Mathematical Society} 33, 351-367.

  Berthelsen, K.K. and \ifelse{latex}{\out{M\o ller}}{Moller}, J. (2003)
  Likelihood and non-parametric Bayesian MCMC inference 
  for spatial point processes based on perfect simulation and
  path sampling. 
  \emph{Scandinavian Journal of Statistics} 30, 549-564.

  \ifelse{latex}{\out{M\o ller}}{Moller}, J. and Waagepetersen, R. (2003).
  \emph{Statistical Inference and Simulation for Spatial Point Processes.}
  Chapman and Hall/CRC.
}
\author{
  \adrian
  
  
  based on original code for the Strauss process by 
  Kasper Klitgaard Berthelsen.
}
\examples{
   X <- rHardcore(0.05,1.5,square(141.4))
   Z <- rHardcore(100,0.05)
}
\seealso{
  \code{\link{rmh}},
  \code{\link{Hardcore}},
  \code{\link{rStrauss}},
  \code{\link{rStraussHard}},
  \code{\link{rDiggleGratton}}.
  \code{\link{rDGS}},
  \code{\link{rPenttinen}}.
}
\keyword{spatial}
\keyword{datagen}
