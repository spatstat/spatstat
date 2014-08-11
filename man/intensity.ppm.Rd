\name{intensity.ppm}
\alias{intensity.ppm}
\title{
  Intensity of Fitted Point Process Model
}
\description{
  Computes the intensity of a fitted point process model.
}
\usage{
\method{intensity}{ppm}(X, ...)
}
\arguments{
  \item{X}{
    A fitted point process model (object of class \code{"ppm"}).
  }
  \item{\dots}{
    Arguments passed to \code{\link{predict.ppm}} in some cases.
    See Details.
  }
}
\details{
  This is a method for the generic function \code{\link{intensity}}
  for fitted point process models (class \code{"ppm"}).

  The intensity of a point process model is the expected
  number of random points per unit area. 

  If \code{X} is a Poisson point process model, the intensity of the
  process is computed exactly.
  The result is a numerical value if \code{X}
  is a stationary Poisson point process, and a pixel image if \code{X}
  is non-stationary. (In the latter case, the resolution of the pixel
  image is controlled by the arguments \code{\dots} which are passed
  to \code{\link{predict.ppm}}.)

  If \code{X} is another Gibbs point process model, the intensity is
  computed approximately using the Poisson-saddlepoint approximation
  (Baddeley and Nair, 2012). Currently this is implemented only for
  pairwise interactions.
}
\value{
  A numeric value (if the model is stationary)
  or a pixel image.
}
\references{
  Baddeley, A. and Nair, G. (2012) 
  Fast approximation of the intensity of Gibbs point processes.
  \emph{Electronic Journal of Statistics} \bold{6} 1155--1169.
}
\seealso{
  \code{\link{intensity}},
  \code{\link{intensity.ppp}}
}
\examples{
  fitP <- ppm(swedishpines, ~1, Poisson())
  intensity(fitP)
  fitS <- ppm(swedishpines, ~1, Strauss(9))
  intensity(fitS)
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Gopal Nair.
}
\keyword{spatial}
\keyword{models}
