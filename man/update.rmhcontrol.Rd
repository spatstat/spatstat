\name{update.rmhcontrol}
\alias{update.rmhcontrol}
\title{Update Control Parameters of Metropolis-Hastings Algorithm}
\description{
  \code{update} method for class \code{"rmhcontrol"}.
}
\usage{
  \method{update}{rmhcontrol}(object, \dots)
}
\arguments{
  \item{object}{
    Object of class \code{"rmhcontrol"} containing control parameters
    for a Metropolis-Hastings algorithm.
  }
  \item{\dots}{
    Arguments to be updated in the new call to \code{\link{rmhcontrol}}.
  }
}
\details{
  This is a method for the generic function \code{\link{update}}
  for the class \code{"rmhcontrol"}. An object of class \code{"rmhcontrol"}
  describes a set of control parameters for the Metropolis-Hastings
  simulation algorithm. See \code{\link{rmhcontrol}}).

  \code{update.rmhcontrol} will modify the parameters
  specified by \code{object} according to the new arguments given.
}
\value{
  Another object of class \code{"rmhcontrol"}.
}
\examples{
  a <- rmhcontrol(expand=1)
  update(a, expand=2)
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{methods}
\keyword{models}
