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

  \method{pseudoR2}{ppm}(object, \dots, keepoffset=TRUE)

  \method{pseudoR2}{lppm}(object, \dots, keepoffset=TRUE)
}
\arguments{
  \item{object}{
    Fitted point process model. An object of class \code{"ppm"}
    or \code{"lppm"}.
  }
  \item{keepoffset}{
    Logical value indicating whether to retain offset terms in the model
    when computing the deviance difference. See Details.
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
  and \eqn{D_0}{D0} is the deviance of the null model.
  Deviance is defined as twice the negative log-likelihood
  or log-pseudolikelihood.

  The null model is usually obtained by re-fitting the model
  using the trend formula \code{~1}.
  However if the original model formula included \code{offset} terms,
  and if \code{keepoffset=TRUE} (the default),
  then the null model formula consists of these offset terms. This
  ensures that the \code{pseudoR2} value is non-negative.
}
\value{
  A single numeric value.
}
\author{
  \spatstatAuthors.
}
\seealso{
  \code{\link{deviance.ppm}}, 
  \code{\link{deviance.lppm}}.
}
\examples{
  fit <- ppm(swedishpines ~ x+y)
  pseudoR2(fit)

  xcoord <- as.im(function(x,y) x, Window(swedishpines))
  fut <- ppm(swedishpines ~ offset(xcoord/200) + y)
  pseudoR2(fut)
}
\keyword{spatial}
\keyword{models}

