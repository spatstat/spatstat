\name{print.ppm}
\alias{print.ppm}
\title{Print a Fitted Point Process Model}
\description{
  Default \code{print} method for a fitted point process model.
}
\usage{
 \method{print}{ppm}(x,\dots,
     what=c("all", "model", "trend", "interaction", "se", "errors"))
}
\arguments{
  \item{x}{
    A fitted point process model, typically obtained from
    the model-fittingg algorithm \code{\link{ppm}}.
    An object of class \code{"ppm"}.
  }
  \item{what}{
    Character vector (partially-matched) indicating what information
    should be printed.
  }
  \item{\dots}{Ignored.}
}
\value{
  none.
}
\details{
  This is the \code{print} method for the class \code{"ppm"}.
  It prints information about the fitted model in a sensible format.

  The argument \code{what} makes it possible to print only some
  of the information.

  If \code{what} is missing, then by default,
  standard errors for the estimated coefficients of the model
  will be printed only if the model is a Poisson point process.
  To print the standard errors for a non-Poisson model,
  call \code{print.ppm} with the argument \code{what} given explicitly,
  or reset the default rule by typing
  \code{spatstat.options(print.ppm.SE="always")}.
}
\seealso{
  \code{\link{ppm.object}} for details of the class \code{"ppm"}.

  \code{\link{ppm}} for generating these objects.
  
  \code{\link{plot.ppm}},
  \code{\link{predict.ppm}}
}
\examples{
 \dontrun{
 m <- ppm(cells, ~1, Strauss(0.05))
 m
 }
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{print}
\keyword{models}

