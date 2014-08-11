\name{rStraussHard}
\alias{rStraussHard}
\title{Perfect Simulation of the Strauss-Hardcore Process}
\description{
  Generate a random pattern of points, a simulated realisation
  of the Strauss-Hardcore process, using a perfect simulation algorithm.
}
\usage{
  rStraussHard(beta, gamma = 1, R = 0, H = 0, W = owin())
}
\arguments{
  \item{beta}{
    intensity parameter (a positive number).
  }
  \item{gamma}{
    interaction parameter (a number between 0 and 1, inclusive).
  }
  \item{R}{
    interaction radius (a non-negative number).
  }
  \item{H}{
    hard core distance (a non-negative number smaller than \code{R}).
  }
  \item{W}{
    window (object of class \code{"owin"}) in which to
    generate the random pattern. Currently this must be a rectangular
    window.
  }
}
\details{
  This function generates a realisation of the 
  Strauss-Hardcore point process in the window \code{W}
  using a \sQuote{perfect simulation} algorithm.

  The Strauss-Hardcore process is described in \code{\link{StraussHard}}.

  The simulation algorithm used to generate the point pattern
  is \sQuote{dominated coupling from the past}
  as implemented by Berthelsen and Moller (2002, 2003).
  This is a \sQuote{perfect simulation} or \sQuote{exact simulation}
  algorithm, so called because the output of the algorithm is guaranteed
  to have the correct probability distribution exactly (unlike the
  Metropolis-Hastings algorithm used in \code{\link{rmh}}, whose output
  is only approximately correct).

  A limitation of the perfect simulation algorithm
  is that the interaction parameter
  \eqn{\gamma}{gamma} must be less than or equal to \eqn{1}.
  To simulate a Strauss-hardcore process with
  \eqn{\gamma > 1}{gamma > 1}, use \code{\link{rmh}}.

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

  Moller, J. and Waagepetersen, R. (2003).
  \emph{Statistical Inference and Simulation for Spatial Point Processes.}
  Chapman and Hall/CRC.
}
\author{
  Kasper Klitgaard Berthelsen and Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
}
\examples{
   Z <- rStraussHard(100,0.7,0.05,0.02)
}
\seealso{
  \code{\link{rmh}},
  \code{\link{rStrauss}},
  \code{\link{StraussHard}}.
}

\keyword{spatial}
\keyword{datagen}
