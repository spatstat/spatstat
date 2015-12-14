\name{simulate.lppm}
\alias{simulate.lppm}
\title{Simulate a Fitted Point Process Model on a Linear Network}
\description{
  Generates simulated realisations from a fitted Poisson
  point process model on a linear network.
}
\usage{
  \method{simulate}{lppm}(object, nsim=1, ...,
                         new.coef=NULL,
                         progress=(nsim > 1),
                         drop=FALSE)
}
\arguments{
  \item{object}{
    Fitted point process model on a linear network.
    An object of class \code{"lppm"}.
  }
  \item{nsim}{
    Number of simulated realisations.
  }
  \item{progress}{
    Logical flag indicating whether to print progress reports
    for the sequence of simulations.
  }
  \item{new.coef}{
    New values for the canonical parameters of the model.
    A numeric vector of the same length as \code{coef(object)}.
  }
  \item{\dots}{
    Arguments passed to \code{\link{predict.lppm}}
    to determine the spatial resolution of the image of the fitted intensity
    used in the simulation.
  }
  \item{drop}{
    Logical. If \code{nsim=1} and \code{drop=TRUE}, the
    result will be a point pattern, rather than a list 
    containing a point pattern.
  }
}
\details{
  This function is a method for the generic function
  \code{\link[stats]{simulate}} for the class \code{"lppm"} of fitted
  point process models on a linear network.

  Only Poisson process models are supported so far.
  
  Simulations are performed by \code{\link{rpoislpp}}.
}
\value{
  A list of length \code{nsim} containing simulated point patterns
  (objects of class \code{"lpp"}) on the same linear network as the
  original data used to fit the model.
  The result also belongs to the class \code{"solist"}, so that it can be
  plotted, and the class \code{"timed"}, so that the total computation
  time is recorded.
}
\examples{
  fit <- lppm(unmark(chicago) ~ y)
  simulate(fit)[[1]]
}
\seealso{
  \code{\link{lppm}},
  \code{\link{rpoislpp}},
  \code{\link[stats]{simulate}}
}
\author{Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  ,
  Rolf Turner \email{r.turner@auckland.ac.nz}
  
  and Ege Rubak \email{rubak@math.aau.dk}
}
\keyword{spatial}
\keyword{models}
