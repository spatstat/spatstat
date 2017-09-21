\name{pseudoR2}
\alias{pseudoR2}
\alias{pseudoR2.ppm}
\alias{pseudoR2.lppm}
\title{
  Calculate Pseudo-R-Squared for Point Process Model
}
\description{
  Given a fitted point process model, calculate 
  the pseudo-R-squared value, which measures the 
  fraction of variation in the data that is explained
  by the model.  
}
\usage{
  pseudoR2(object, \dots)

  \method{pseudoR2}{ppm}(object, \dots)

  \method{pseudoR2}{lppm}(object, \dots)
}
\arguments{
  \item{object}{
    Fitted point process model. An object of class \code{"ppm"}
    or \code{"lppm"}.
  }
  \item{\dots}{
    Additional arguments passed to \code{\link{deviance.ppm}}
    or \code{\link{deviance.lppm}}.
  }
}
\details{
  The function \code{pseudoR2} is generic, with methods
  for fitted point process models of class \code{"ppm"} and
  \code{"lppm"}.

  This function computes McFadden's pseudo-Rsquared
  \deqn{
    R^2 = 1 - \frac{D}{D_0}
  }{
    R^2 = 1 - D/D0
  }
  where \eqn{D} is the deviance of the fitted model \code{object},
  and \eqn{D_0}{D0} is the deviance of the null model
  (obtained by refitting \code{object}
  using the trend formula \code{~1}).
  Deviance is defined as twice the negative log-likelihood
  or log-pseudolikelihood.
}
\value{
  A single numeric value.
}
\author{
  \adrian
  
  
  \rolf
  
  and \ege
  
}
\seealso{
  \code{\link{deviance.ppm}}, 
  \code{\link{deviance.lppm}}.
}
\examples{
  fit <- ppm(swedishpines ~ x+y)
  pseudoR2(fit)
}
\keyword{spatial}
\keyword{models}

