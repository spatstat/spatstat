\name{parameters}
\alias{parameters}
\alias{parameters.kppm}
\alias{parameters.ppm}
\alias{parameters.profilepl}
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

\method{parameters}{kppm}(model, \dots)

\method{parameters}{ppm}(model, \dots)

\method{parameters}{profilepl}(model, \dots)
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
  for class \code{"ppm"}, \code{"kppm"} and \code{"profilepl"}
  and possibly for other classes.
}
\value{
  A named list, whose format depends on the fitted model.
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
