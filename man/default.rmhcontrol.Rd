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
   default.rmhcontrol(model)
}
\arguments{
  \item{model}{
    A fitted point process model (object of class \code{"ppm"})
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
  \code{\link{default.expand}(model)}.

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
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{datagen}
