\name{parameters}
\alias{parameters}
\alias{parameters.dppm}
\alias{parameters.kppm}
\alias{parameters.ppm}
\alias{parameters.profilepl}
\alias{parameters.interact}
\alias{parameters.fii}
\title{
  Extract Model Parameters in Understandable Form
}
\description{
  Given a fitted model of some kind, this function
  extracts all the parameters needed to specify the model,
  and returns them as a list.
}
\usage{
parameters(model, \dots)

\method{parameters}{dppm}(model, \dots)

\method{parameters}{kppm}(model, \dots)

\method{parameters}{ppm}(model, \dots)

\method{parameters}{profilepl}(model, \dots)

\method{parameters}{fii}(model, \dots)

\method{parameters}{interact}(model, \dots)
}
\arguments{
  \item{model}{
    A fitted model of some kind.
  }
  \item{\dots}{
    Arguments passed to methods.
  }
}
\details{
  The argument \code{model} should be a fitted model of some kind.
  This function extracts all the parameters that would be needed to
  specify the model, and returns them as a list.

  The function \code{parameters} is generic, with methods
  for class \code{"ppm"}, \code{"kppm"}, \code{"dppm"} and \code{"profilepl"}
  and other classes.
}
\value{
  A named list, whose format depends on the fitted model.
}
\author{
\spatstatAuthors
}
\seealso{
  \code{\link{coef}}
}
\examples{
  fit1 <- ppm(cells ~ x, Strauss(0.1))
  parameters(fit1)
  fit2 <- kppm(redwood ~ x, "Thomas")
  parameters(fit2)
}
\keyword{spatial}
\keyword{models}
