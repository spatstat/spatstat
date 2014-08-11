\name{as.ppm}
\alias{as.ppm}
\alias{as.ppm.ppm}
\alias{as.ppm.profilepl}
\alias{as.ppm.kppm}
\title{Extract Fitted Point Process Model}
\description{
  Extracts the fitted point process model from
  some kind of fitted model.
}
\usage{
as.ppm(object)

\method{as.ppm}{ppm}(object)

\method{as.ppm}{profilepl}(object)

\method{as.ppm}{kppm}(object)
}
\arguments{
  \item{object}{An object that includes a
    fitted Poisson or Gibbs point process model.
    An object of class \code{"ppm"}, \code{"profilepl"} or \code{"kppm"}
    or possibly other classes.
  }
}
\details{
  The function \code{as.ppm} extracts 
  the fitted point process model (of class \code{"ppm"})
  from a suitable object.
  
  The function \code{as.ppm} is generic, with methods for the classes
  \code{"ppm"}, \code{"profilepl"} and \code{"kppm"},
  and possibly for other classes.

  For the class \code{"profilepl"} of models fitted by maximum profile
  pseudolikelihood, the method \code{as.ppm.profilepl} extracts the
  fitted point process model (with the optimal values of the
  irregular parameters).

  For the class \code{"kppm"} of models fitted by minimum contrast
  using Waagepetersen's two-step estimation procedure
  (see \code{\link{kppm}}), the method \code{as.ppm.kppm}
  extracts the Poisson point process model that is fitted in the
  first stage of the procedure.
}
\value{
  An object of class \code{"ppm"}.
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\seealso{
  \code{\link{ppm}},
  \code{\link{profilepl}}.
}
\examples{
   # fit a model by profile maximum pseudolikelihood
   rvals <- data.frame(r=(1:10)/100)
   pfit <- profilepl(rvals, Strauss, cells, ~1)
   # extract the fitted model
   fit <- as.ppm(pfit)
}
\keyword{spatial}
\keyword{models}
