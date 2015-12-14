\name{default.expand}
\alias{default.expand}
\title{Default Expansion Rule for Simulation of Model}
\description{
  Defines the default expansion window or expansion rule
  for simulation of a fitted point process model.
}
\usage{
  default.expand(object, m=2, epsilon=1e-6, w=Window(object))
}
\arguments{
  \item{object}{
    A point process model (object of class \code{"ppm"}
    or \code{"rmhmodel"}).
  }
  \item{m}{
    A single numeric value.
    The window will be expanded by a distance 
    \code{m * reach(object)} along each side.
  }
  \item{epsilon}{
    Threshold argument passed to \code{\link{reach}} to determine
    \code{reach(object)}.
  }
  \item{w}{
    Optional. The un-expanded window in which the model is defined.
    The resulting simulated point patterns will lie in this window. 
  }
} 
\value{
  A window expansion rule (object of class \code{"rmhexpand"}).
}
\details{
  This function computes a default value for the
  expansion rule (the argument \code{expand} in \code{\link{rmhcontrol}})
  given a fitted point process model \code{object}.
  This default is used by \code{\link{envelope}}, 
  \code{\link{qqplot.ppm}}, \code{\link{simulate.ppm}} and other functions.

  Suppose we wish to generate simulated realisations
  of a fitted point process model inside a window \code{w}.
  It is advisable to first simulate
  the pattern on a larger window, and then clip it to the original
  window \code{w}. This avoids edge effects in the simulation.
  It is called \emph{expansion} of the simulation window.
  
  Accordingly, for the Metropolis-Hastings simulation algorithm
  \code{\link{rmh}}, the algorithm control parameters specified by
  \code{\link{rmhcontrol}} include an argument \code{expand} that
  determines the expansion of the simulation window.

  The function \code{default.expand} determines the default expansion
  rule for a fitted point process model \code{object}.
  
  If the model is Poisson, then no expansion is necessary.
  No expansion is performed by default,
  and \code{default.expand} returns a rule representing no expansion.
  The simulation window is the original window \code{w = Window(object)}.
  
  If the model depends on external covariates (i.e.\ covariates other than
  the Cartesian covariates \code{x} and \code{y} and the \code{marks})
  then no expansion is feasible, in general, because the spatial domain
  of the covariates is not guaranteed to be large enough.
  \code{default.expand} returns a rule representing no expansion.
  The simulation window is the original window \code{w = Window(object)}.

  If the model depends on the Cartesian covariates \code{x} and \code{y},
  it would be feasible to expand the simulation window, and this was the
  default for \pkg{spatstat} version 1.24-1 and earlier.
  However this sometimes produces artefacts (such as an empty point pattern)
  or memory overflow, because the fitted trend, extrapolated outside the
  original window of the data, may become very large.
  In \pkg{spatstat} version 1.24-2 and later, the
  default rule is \emph{not} to expand if the model depends
  on \code{x} or \code{y}.
  Again \code{default.expand} returns a rule representing no expansion.
  
  Otherwise, expansion will occur.
  The original window \code{w = Window(object)} is expanded by
  a distance \code{m * rr}, where
  \code{rr} is the interaction range of the model, computed by
  \code{\link{reach}}. If \code{w} is a rectangle then
  each edge of \code{w} is displaced outward by distance \code{m * rr}.
  If \code{w} is not a rectangle then \code{w} is dilated by
  distance \code{m * rr} using \code{\link{dilation}}.
}
\seealso{
  \code{\link{rmhexpand}},
  \code{\link{rmhcontrol}},
  \code{\link{rmh}},
  \code{\link{envelope}},
  \code{\link{qqplot.ppm}}
}
\examples{
  data(cells)
  fit <- ppm(cells, ~1, Strauss(0.07))
  default.expand(fit)
  mod <- rmhmodel(cif="strauss", par=list(beta=100, gamma=0.5, r=0.07))
  default.expand(fit)
}
\author{Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  and Rolf Turner \email{r.turner@auckland.ac.nz}
  
}
\keyword{spatial}
\keyword{datagen}
