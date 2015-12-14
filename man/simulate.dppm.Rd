\name{simulate.dppm}
\alias{simulate.dppm}
\alias{simulate.detpointprocfamily}
\title{Simulation of Determinantal Point Process Model}
\description{
  Generates simulated realisations from a determinantal point process model.
}
\usage{
  \method{simulate}{dppm}(object, nsim = 1, seed = NULL, \dots,
    W = NULL, trunc = 0.99, correction = "periodic", rbord = reach(object))

  \method{simulate}{detpointprocfamily}(object, nsim = 1, seed = NULL, \dots,
    W = NULL, trunc = 0.99, correction = "periodic", rbord = reach(object))
}
\arguments{
  \item{object}{
    Determinantal point process model. An object of class
    \code{"detpointprocfamily"} or \code{"dppm"}.
  }
  \item{nsim}{Number of simulated realisations.}
  \item{seed}{
    an object specifying whether and how to initialise the random
    number generator. Either \code{NULL} or an integer that will be
    used in a call to \code{\link[base:Random]{set.seed}} before
    simulating the point patterns.
  }
  \item{\dots}{Arguments passed on to \code{\link{rdpp}}.}
  \item{W}{
    Object specifying the window of simulation (defaults to a unit
    box if nothing else is sensible -- see Details). Can be any single
    argument acceptable to \code{\link{as.boxx}} (e.g. an \code{"owin"},
    \code{"box3"} or \code{"boxx"} object).
  }
  \item{trunc}{
    Numeric value specifying how the model truncation is preformed. See
    Details.
  }
  \item{correction}{
    Character string specifying the type of correction to use.
    The options are "periodic" (default) and "border". See Details.
  }
  \item{rbord}{
    Numeric value specifying the extent of the border correction if this
    correction is used. See Details.
  }
}
\details{
  These functions are methods for the generic function
  \code{\link[stats]{simulate}} for the classes \code{"detpointprocfamily"} and
  \code{"dppm"} of determinantal point process models.

  The return value is a list of \code{nsim} point patterns.
  It also carries an attribute \code{"seed"} that
  captures the initial state of the random number generator.
  This follows the convention used in
  \code{simulate.lm} (see \code{\link[stats]{simulate}}).
  It can be used to force a sequence of simulations to be
  repeated exactly, as shown in the examples for
  \code{\link[stats]{simulate}}.

  The exact simulation of a determinantal point process model involves
  an infinite series, which typically has no analytical solution. In the
  implementation a truncation is performed. The truncation
  \code{trunc} can be specified either directly as a positive integer
  or as a fraction between 0 and 1. In the latter case the truncation is chosen
  such that the expected number of points in a simulation is
  \code{trunc} times the theoretical expected number of points in the
  model. The default is 0.99.

  The window of the returned point pattern(s) can be specified via the
  argument \code{W}. For a fitted model (of class \code{"dppm"}) it
  defaults to the observation window of the data used to fit the
  model. For inhomogeneous models it defaults to the window of the
  intensity image. Otherwise it defaults to a unit box.  For
  non-rectangular windows simulation is done in the containing rectangle
  and then restricted to the window.  For inhomogeneous models a
  stationary model is first simulated using the maximum intensity and
  then the result is obtained by thinning.

  The default is to use periodic edge correction for simulation such
  that opposite edges are glued together.  If border correction is used
  then the simulation is done in an extended window. Edge effects are
  theoretically completely removed by doubling the size of the window in
  each spatial dimension, but for practical purposes much less extension
  may be sufficient. The numeric \code{rbord} determines the extend of
  the extra space added to the window.
}
\value{
  A list of length \code{nsim} containing simulated point patterns
  (objects of class \code{"ppp"}). The list has class \code{"solist"}.

  The return value also carries an attribute \code{"seed"} that
  captures the initial state of the random number generator.
  See Details.
}
\references{
  Lavancier, F. \ifelse{latex}{\out{M\o ller}}{Moller}, J. and Rubak, E. (2015)
  Determinantal point process models and statistical inference
  \emph{Journal of the Royal Statistical Society, Series B}
  \bold{77}, 853--977.
}
\author{
  Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  Rolf Turner \email{r.turner@auckland.ac.nz}
  
  and Ege Rubak \email{rubak@math.aau.dk}
  
}
\seealso{
  \code{\link{rdpp}},
  \code{\link[stats]{simulate}}
}
\examples{
model <- dppGauss(lambda=100, alpha=.05, d=2)
simulate(model, 2)
}
\keyword{datagen}
\keyword{spatial}
\keyword{models}
