\name{rStrauss}
\alias{rStrauss}
\title{Perfect Simulation of the Strauss Process}
\description{
  Generate a random pattern of points, a simulated realisation
  of the Strauss process, using a perfect simulation algorithm.
}
\usage{
  rStrauss(beta, gamma = 1, R = 0, W = owin(), expand=TRUE, nsim=1, drop=TRUE)
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
  Strauss point process in the window \code{W}
  using a \sQuote{perfect simulation} algorithm.

  The Strauss process (Strauss, 1975; Kelly and Ripley, 1976)
  is a model for spatial inhibition, ranging from
  a strong `hard core' inhibition to a completely random pattern
  according to the value of \code{gamma}.

  The Strauss process with interaction radius \eqn{R} and 
  parameters \eqn{\beta}{beta} and \eqn{\gamma}{gamma}
  is the pairwise interaction point process
  with probability density 
  \deqn{
    f(x_1,\ldots,x_n) =
    \alpha \beta^{n(x)} \gamma^{s(x)}
  }{
    f(x_1,\ldots,x_n) =
    alpha . beta^n(x) gamma^s(x)
  }
  where \eqn{x_1,\ldots,x_n}{x[1],\ldots,x[n]} represent the 
  points of the pattern, \eqn{n(x)} is the number of points in the
  pattern, \eqn{s(x)} is the number of distinct unordered pairs of
  points that are closer than \eqn{R} units apart,
  and \eqn{\alpha}{alpha} is the normalising constant.
  Intuitively, each point of the pattern
  contributes a factor \eqn{\beta}{beta} to the 
  probability density, and each pair of points
  closer than \eqn{r} units apart contributes a factor
  \eqn{\gamma}{gamma} to the density.

  The interaction parameter \eqn{\gamma}{gamma} must be less than
  or equal to \eqn{1} in order that the process be well-defined
  (Kelly and Ripley, 1976).
  This model describes an ``ordered'' or ``inhibitive'' pattern.
  If \eqn{\gamma=1}{gamma=1} it reduces to a Poisson process
  (complete spatial randomness) with intensity \eqn{\beta}{beta}.
  If \eqn{\gamma=0}{gamma=0} it is called a ``hard core process''
  with hard core radius \eqn{R/2}, since no pair of points is permitted
  to lie closer than \eqn{R} units apart.

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

  Kelly, F.P. and Ripley, B.D. (1976)
  On Strauss's model for clustering.
  \emph{Biometrika} \bold{63}, 357--360.

  \ifelse{latex}{\out{M\o ller}}{Moller}, J. and Waagepetersen, R. (2003).
  \emph{Statistical Inference and Simulation for Spatial Point Processes.}
  Chapman and Hall/CRC.

  Strauss, D.J. (1975)
  A model for clustering.
  \emph{Biometrika} \bold{62}, 467--475.
}
\author{
  Kasper Klitgaard Berthelsen,
  adapted for \pkg{spatstat} by \adrian
}
\examples{
   X <- rStrauss(0.05,0.2,1.5,square(141.4))
   Z <- rStrauss(100,0.7,0.05)
}
\seealso{
  \code{\link{rmh}},
  \code{\link{Strauss}},
  \code{\link{rHardcore}},
  \code{\link{rStraussHard}},
  \code{\link{rDiggleGratton}},
  \code{\link{rDGS}},
  \code{\link{rPenttinen}}.
}

\keyword{spatial}
\keyword{datagen}
