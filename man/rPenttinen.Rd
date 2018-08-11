\name{rPenttinen}
\alias{rPenttinen}
\title{Perfect Simulation of the Penttinen Process}
\description{
  Generate a random pattern of points, a simulated realisation
  of the Penttinen process, using a perfect simulation algorithm.
}
\usage{
  rPenttinen(beta, gamma=1, R, W = owin(), expand=TRUE, nsim=1, drop=TRUE)
}
\arguments{
  \item{beta}{
    intensity parameter (a positive number).
  }
  \item{gamma}{
    Interaction strength parameter (a number between 0 and 1).
  }
  \item{R}{
    disc radius (a non-negative number).
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
  Penttinen point process in the window \code{W}
  using a \sQuote{perfect simulation} algorithm.

  Penttinen (1984, Example 2.1, page 18), citing Cormack (1979),
  described the pairwise interaction point process with interaction factor
  \deqn{
    h(d) = e^{\theta A(d)} = \gamma^{A(d)}
  }{
    h(d) = exp(theta * A(d)) = gamma^(A(d))
  }
  between each pair of points separated by a distance $d$.
  Here \eqn{A(d)} is the area of intersection between two discs
  of radius \eqn{R} separated by a distance \eqn{d}, normalised so that
  \eqn{A(0) = 1}.

  The simulation algorithm used to generate the point pattern
  is \sQuote{dominated coupling from the past}
  as implemented by
  Berthelsen and \ifelse{latex}{\out{M\o ller}}{Moller} (2002, 2003).
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

  Cormack, R.M. (1979)
  Spatial aspects of competition between individuals.
  Pages 151--212 in \emph{Spatial and Temporal Analysis in Ecology},
  eds. R.M. Cormack and J.K. Ord, International Co-operative
  Publishing House, Fairland, MD, USA. 

  \ifelse{latex}{\out{M\o ller}}{Moller}, J. and Waagepetersen, R. (2003).
  \emph{Statistical Inference and Simulation for Spatial Point Processes.}
  Chapman and Hall/CRC.

  Penttinen, A. (1984) 
  \emph{Modelling Interaction in Spatial Point Patterns:
  Parameter Estimation by the Maximum Likelihood Method.}
  \ifelse{latex}{\out{Jyv\"askyl\"a}}{Jyvaskyla}
  Studies in Computer Science, Economics and Statistics \bold{7},
  University of \ifelse{latex}{\out{Jyv\"askyl\"a}}{Jyvaskyla}, Finland.
}
\author{
  \adrian,
  based on original code for the Strauss process by 
  Kasper Klitgaard Berthelsen.
}
\examples{
   X <- rPenttinen(50, 0.5, 0.02)
}
\seealso{
  \code{\link{rmh}},
  \code{\link{Penttinen}}.

  \code{\link{rStrauss}},
  \code{\link{rHardcore}},
  \code{\link{rStraussHard}},
  \code{\link{rDiggleGratton}},
  \code{\link{rDGS}}.
}
\keyword{spatial}
\keyword{datagen}
