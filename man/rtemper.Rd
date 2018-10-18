\name{rtemper}
\alias{rtemper}
\title{
  Simulated Annealing or Simulated Tempering for Gibbs Point Processes
}
\description{
  Performs simulated annealing or simulated tempering
  for a Gibbs point process model using a specified
  annealing schedule.
}
\usage{
rtemper(model, invtemp, nrep, \dots, track=FALSE, start = NULL, verbose = FALSE)
}
\arguments{
  \item{model}{
    A Gibbs point process model: a fitted Gibbs point process model
    (object of class \code{"ppm"}), or any data acceptable to
    \code{\link{rmhmodel}}.
  }
  \item{invtemp}{
    A numeric vector of positive numbers.
    The sequence of values of inverse temperature that will be used.
  }
  \item{nrep}{
    An integer vector of the same length as \code{invtemp}.
    The value \code{nrep[i]} specifies the number of steps of the
    Metropolis-Hastings algorithm that will be performed
    at inverse temperature \code{invtemp[i]}.
  }
  \item{start}{
    Initial starting state for the simulation.
    Any data acceptable to \code{\link{rmhstart}}.
  }
  \item{track}{
    Logical flag indicating whether to save the transition
    history of the simulations.
  }
  \item{\dots}{
    Additional arguments passed to \code{\link{rmh.default}}.
  }
  \item{verbose}{
    Logical value indicating whether to print progress reports.
  }
}
\details{
  The Metropolis-Hastings simulation algorithm
  \code{\link{rmh}} is run for
  \code{nrep[1]} steps at inverse temperature \code{invtemp[1]},
  then for 
  \code{nrep[2]} steps at inverse temperature \code{invtemp[2]},
  and so on.

  Setting the inverse temperature to a value \eqn{\alpha}{alpha}
  means that the probability density of the Gibbs model, \eqn{f(x)},
  is replaced by \eqn{g(x) = C\, f(x)^\alpha}{g(x) = C f(x)^alpha}
  where \eqn{C} is a normalising constant depending on
  \eqn{\alpha}{alpha}.
  Larger values of \eqn{\alpha}{alpha} exaggerate the high and low
  values of probability density, while smaller values of \eqn{\alpha}{alpha}
  flatten out the probability density.

  For example if the original \code{model} is a Strauss process,
  the modified model is close to a hard core process
  for large values of inverse temperature, and close to a Poisson process
  for small values of inverse temperature. 
}
\value{
  A point pattern (object of class \code{"ppp"}).

  If \code{track=TRUE}, the result also has an attribute
  \code{"history"} which is a data frame with columns
  \code{proposaltype}, \code{accepted}, 
  \code{numerator} and \code{denominator}, as described
  in \code{\link{rmh.default}}.
}
\author{
  \spatstatAuthors.
}
\seealso{
  \code{\link{rmh.default}},
  \code{\link{rmh}}.
}
\examples{
   stra <- rmhmodel(cif="strauss",
                    par=list(beta=2,gamma=0.2,r=0.7),
                    w=square(10))
   nr <- if(interactive()) 1e5 else 1e4
   Y <- rtemper(stra, c(1, 2, 4, 8), nr * (1:4), verbose=TRUE, track=TRUE)
}
\keyword{spatial}
\keyword{datagen}
