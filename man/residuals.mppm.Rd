\name{residuals.mppm}
\alias{residuals.mppm}
\title{Residuals for Point Process Model Fitted to Multiple Point Patterns}
\description{
  Given a point process model fitted to multiple point patterns,
  compute residuals for each pattern.
}
\usage{
\method{residuals}{mppm}(object, type = "raw", ..., 
                             fittedvalues = fitted.mppm(object))
}
\arguments{
  \item{object}{Fitted point process model (object of class \code{"mppm"}).}
  \item{\dots}{Ignored.}
  \item{type}{Type of residuals: either \code{"raw"}, \code{"pearson"}
    or \code{"inverse"}. Partially matched.}
  \item{fittedvalues}{Advanced use only.
    Fitted values of the model to be used in the calculation.
  }
}
\details{
  Baddeley et al (2005) defined residuals for the fit of
  a point process model to spatial point pattern data.
  For an explanation of these residuals, see the help file for
  \code{\link[spatstat]{residuals.ppm}}.

  This function computes the residuals
  for a point process model fitted to \emph{multiple} point patterns.
  The \code{object} should be an object of class \code{"mppm"} obtained
  from \code{\link{mppm}}.

  The return value is a list.
  The number of entries in the list equals the
  number of point patterns in the original data. Each entry in the list
  has the same format as the output of
  \code{\link[spatstat]{residuals.ppm}}.
  That is, each entry in the list is a signed measure (object of
  class \code{"msr"}) giving the residual measure for the corresponding
  point pattern. 
}
\value{
  A list of signed measures (objects of class \code{"msr"})
  giving the residual measure for each of the original point patterns.
  See Details.
}
\examples{
    fit <- mppm(Bugs ~ x, hyperframe(Bugs=waterstriders))
    r <- residuals(fit)
    # compute total residual for each point pattern
    rtot <- sapply(r, integral.msr)
    # standardise the total residuals
    areas <- sapply(windows.mppm(fit), area.owin)
    rtot/sqrt(areas)
}
\references{
  Baddeley, A., Turner, R., Moller, J. and Hazelton, M. (2005)
  Residual analysis for spatial point processes.
  \emph{Journal of the Royal Statistical Society, Series B}
  \bold{67}, 617--666.
}
\author{Adrian Baddeley
  \email{adrian.baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{rolf@math.unb.ca}
  \url{http://www.math.unb.ca/~rolf}
}
\seealso{
  \code{\link{mppm}},
  \code{\link{residuals.mppm}},
  \code{\link{mppm}}
}
\keyword{spatial}
\keyword{models}
