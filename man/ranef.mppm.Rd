\name{ranef.mppm}
\alias{ranef.mppm}
\title{
  Extract Random Effects from Point Process Model
}
\description{
  Given a point process model fitted to a list of point patterns,
  extract the fixed effects of the model.
  A method for \code{ranef}.
}
\usage{
 \method{ranef}{mppm}(object, \dots)
}
\arguments{
  \item{object}{
    A fitted point process model (an object of class \code{"mppm"}).
  }
  \item{\dots}{
    Ignored.
  }
}
\details{
  This is a method for the generic function \code{\link[nlme]{ranef}}.

  The argument \code{object} must be a fitted point process model
  (object of class \code{"mppm"}) produced by the 
  fitting algorithm \code{\link{mppm}}). This represents a
  point process model that has been fitted
  to a list of several point pattern datasets. See \code{\link{mppm}}
  for information.

  This function extracts the coefficients of the random effects
  of the model.
}
\value{
  A data frame, or list of data frames, as described in the help for
  \code{\link[nlme]{ranef.lme}}.
}
\references{
  Baddeley, A., Rubak, E. and Turner, R. (2015)
  \emph{Spatial Point Patterns: Methodology and Applications with R}.
  London: Chapman and Hall/CRC Press. 
}
\author{
  Adrian Baddeley, Ida-Maria Sintorn and Leanne Bischoff.
  Implemented in \pkg{spatstat} by
  \spatstatAuthors.
}
\seealso{
  \code{\link{fixef.mppm}},
  \code{\link{coef.mppm}}
}
\examples{
 H <- hyperframe(Y = waterstriders)
 # Tweak data to exaggerate differences
 H$Y[[1]] <- rthin(H$Y[[1]], 0.3)

 m1 <- mppm(Y ~ id,  data=H, Strauss(7))
 ranef(m1)
 m2 <- mppm(Y ~ 1,  random=~1|id, data=H, Strauss(7))
 ranef(m2)
}
\keyword{spatial}
\keyword{methods}
\keyword{models}
