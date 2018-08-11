\name{simulate.slrm}
\alias{simulate.slrm}
\title{Simulate a Fitted Spatial Logistic Regression Model}
\description{
  Generates simulated realisations from a fitted
  spatial logistic regresson model
}
\usage{
  \method{simulate}{slrm}(object, nsim = 1, seed=NULL, ...,
         window=NULL, covariates=NULL, verbose=TRUE, drop=FALSE)
}
\arguments{
  \item{object}{
    Fitted spatial logistic regression model. An object of class \code{"slrm"}.
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
  \item{drop}{
    Logical. If \code{nsim=1} and \code{drop=TRUE}, the
    result will be a point pattern, rather than a list 
    containing a point pattern.
  }
}
\details{
  This function is a method for the generic function
  \code{\link[stats]{simulate}} for the class \code{"slrm"} of fitted
  spatial logistic regression models.
  
  Simulations are performed by \code{\link{rpoispp}}
  after the intensity has been computed by \code{\link{predict.slrm}}.

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
  X <- copper$SouthPoints
  fit <- slrm(X ~ 1)
  simulate(fit, 2)
  fitxy <- slrm(X ~ x+y)
  simulate(fitxy, 2, window=square(2))
}
\seealso{
  \code{\link{slrm}},
  \code{\link{rpoispp}},
  \code{\link{simulate.ppm}},
  \code{\link{simulate.kppm}},
  \code{\link[stats]{simulate}}
}
\author{\adrian
  
  
  and \rolf
  
}
\keyword{spatial}
\keyword{models}
