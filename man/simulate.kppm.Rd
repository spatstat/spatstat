\name{simulate.kppm}
\alias{simulate.kppm}
\title{Simulate a Fitted Cluster Point Process Model}
\description{
  Generates simulated realisations from a fitted cluster point process model.
}
\usage{
  \method{simulate}{kppm}(object, nsim = 1, seed=NULL, ...,
         window=NULL, covariates=NULL, verbose=TRUE, retry=10,
         drop=FALSE)
}
\arguments{
  \item{object}{
    Fitted cluster point process model. An object of class \code{"kppm"}.
  }
  \item{nsim}{
    Number of simulated realisations.
  }
  \item{seed}{
    an object specifying whether and how to initialise
    the random number generator. Either \code{NULL} or an integer that will
    be used in a call to \code{\link[base:Random]{set.seed}}
    before simulating the point patterns. 
  }
  \item{\dots}{Ignored.}
  \item{window}{
    Optional. Window (object of class \code{"owin"}) in which the
    model should be simulated.
  }
  \item{covariates}{
    Optional. A named list containing new values for the covariates in the
    model. 
  }
  \item{verbose}{
    Logical. Whether to print progress reports (when \code{nsim > 1}).
  }
  \item{retry}{
    Number of times to repeat the simulation if it fails
    (e.g. because of insufficient memory).
  }
  \item{drop}{
    Logical. If \code{nsim=1} and \code{drop=TRUE}, the
    result will be a point pattern, rather than a list 
    containing a point pattern.
  }
}
\details{
  This function is a method for the generic function
  \code{\link[stats]{simulate}} for the class \code{"kppm"} of fitted
  cluster point process models.
  
  Simulations are performed by \code{\link{rThomas}},
  \code{\link{rMatClust}} or \code{\link{rLGCP}}
  depending on the model.

  The return value is a list of point patterns.
  It also carries an attribute \code{"seed"} that
  captures the initial state of the random number generator.
  This follows the convention used in
  \code{simulate.lm} (see \code{\link[stats]{simulate}}).
  It can be used to force a sequence of simulations to be
  repeated exactly, as shown in the examples for \code{\link[stats]{simulate}}.
}
\value{
  A list of length \code{nsim} containing simulated point patterns
  (objects of class \code{"ppp"}).

  The return value also carries an attribute \code{"seed"} that
  captures the initial state of the random number generator.
  See Details.
}
\examples{
  fit <- kppm(redwood ~1, "Thomas")
  simulate(fit, 2)
  fitx <- kppm(redwood ~x, "Thomas")
  simulate(fitx, 2)
}
\seealso{
  \code{\link{kppm}},
  \code{\link{rThomas}},
  \code{\link{rMatClust}},
  \code{\link{rLGCP}},
  \code{\link{simulate.ppm}},
  \code{\link[stats]{simulate}}
}
\author{
  \spatstatAuthors
}
\keyword{spatial}
\keyword{models}
