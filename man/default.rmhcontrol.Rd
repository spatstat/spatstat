\name{default.rmhcontrol}
\alias{default.rmhcontrol}
\title{Set Default Control Parameters for Metropolis-Hastings Algorithm.}
\description{
  Given a fitted point process model, this command
  sets appropriate default values of the
  parameters controlling the iterative behaviour
  of the Metropolis-Hastings algorithm.
}
\usage{
   default.rmhcontrol(model, w=NULL)
}
\arguments{
  \item{model}{
    A fitted point process model (object of class \code{"ppm"})
  }
  \item{w}{
    Optional. Window for the resulting simulated patterns.
  }
}
\value{
  An object of class \code{"rmhcontrol"}. See \code{\link{rmhcontrol}}.
}
\details{
  This function sets the values of 
  the parameters controlling the iterative behaviour
  of the Metropolis-Hastings
  simulation algorithm. It uses default values
  that would be appropriate for the fitted point
  process model \code{model}.

  The expansion parameter \code{expand} is set to
  \code{\link{default.expand}(model, w)}.

  All other parameters revert to their defaults given in
  \code{\link{rmhcontrol.default}}.
  
  See \code{\link{rmhcontrol}} for the full list of control parameters.
  To override default parameters, use \code{\link{update.rmhcontrol}}.
}
\seealso{
  \code{\link{rmhcontrol}},
  \code{\link{update.rmhcontrol}},
  \code{\link{ppm}},
  \code{\link{default.expand}}
}
\examples{
  fit <- ppm(cells, ~1, Strauss(0.1))
  default.rmhcontrol(fit)
  default.rmhcontrol(fit, w=square(2))
}
\author{Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  and Rolf Turner \email{r.turner@auckland.ac.nz}
  
}
\keyword{spatial}
\keyword{datagen}
