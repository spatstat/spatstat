\name{residuals.dppm}
\alias{residuals.dppm}
\title{
  Residuals for Fitted Determinantal Point Process Model
}
\description{
  Given a determinantal point process model fitted to a point pattern,
  compute residuals.
}
\usage{
  \method{residuals}{dppm}(object, \dots)
}
\arguments{
  \item{object}{
    The fitted determinatal point process model (an object of class \code{"dppm"})
    for which residuals should be calculated.
  }
  \item{\dots}{
    Arguments passed to \code{\link{residuals.ppm}}.
  }
}
\value{
  An object of class \code{"msr"}
  representing a signed measure or vector-valued measure
  (see \code{\link{msr}}). This object can be plotted.
}
\details{
  This function extracts the intensity component of the model using
  \code{\link{as.ppm}} and then applies \code{\link{residuals.ppm}}
  to compute the residuals.

  Use \code{\link{plot.msr}} to plot the residuals directly.
}
\seealso{
 \code{\link{msr}},
 \code{\link{dppm}}
}
\examples{
   fit <- dppm(swedishpines ~ x, dppGauss())
   rr <- residuals(fit)
}
\author{
Adrian Baddeley
\email{Adrian.Baddeley@uwa.edu.au}
\url{http://www.maths.uwa.edu.au/~adrian/}
Rolf Turner
\email{r.turner@auckland.ac.nz}
and Ege Rubak
\email{rubak@math.aau.dk}
}
\keyword{spatial}
\keyword{models}
\keyword{methods}
