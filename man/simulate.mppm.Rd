\name{simulate.mppm}
\alias{simulate.mppm}
\title{Simulate a Point Process Model Fitted to Several Point Patterns}
\description{
  Generates simulated realisations from a 
  point process model that was fitted to several point patterns.
}
\usage{
  \method{simulate}{mppm}(object, nsim=1, \dots, verbose=TRUE)
}
\arguments{
  \item{object}{
    Point process model fitted to several point patterns.
    An object of class \code{"mppm"}.
  }
  \item{nsim}{
    Number of simulated realisations (of each original pattern).
  }
  \item{\dots}{
    Further arguments passed to \code{\link{simulate.ppm}}
    to control the simulation.
  }
  \item{verbose}{
    Logical value indicating whether to print progress reports.
  }
}
\details{
  This function is a method for the generic function
  \code{\link[stats]{simulate}} for the class \code{"mppm"} of fitted
  point process models for replicated point pattern data.

  The result is a hyperframe with \code{n} rows and \code{nsim} columns,
  where \code{n} is the number of original point pattern
  datasets to which the model was fitted. Each column of the hyperframe
  contains a simulated version of the original data.
  
  For each of the original point pattern datasets, the
  fitted model for this dataset is extracted using
  \code{\link{subfits}}, then \code{nsim} simulated realisations
  of this model are generated using \code{\link{simulate.ppm}},
  and these are stored in the corresponding row of the output.
}
\value{
  A hyperframe.
}
\examples{
  H <- hyperframe(Bugs=waterstriders)
  fit <- mppm(Bugs ~ id, H)
  y <- simulate(fit, nsim=2)
  y
  plot(y[1,,drop=TRUE], main="Simulations for Waterstriders pattern 1")
  plot(y[,1,drop=TRUE], main="Simulation 1 for each Waterstriders pattern")
}
\seealso{
  \code{\link{mppm}},
  \code{\link{simulate.ppm}}.
}
\author{
  Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au},
 Rolf Turner \email{r.turner@auckland.ac.nz}
 and Ege Rubak \email{rubak@math.aau.dk}.
}
\keyword{spatial}
\keyword{models}
