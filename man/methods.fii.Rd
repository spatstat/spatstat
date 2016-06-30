\name{methods.fii}
\alias{methods.fii} %DoNotExport
\Rdversion{1.1}
\alias{print.fii}
\alias{plot.fii}
\alias{coef.fii}
\alias{summary.fii}
\alias{print.summary.fii}
\alias{coef.summary.fii}
\title{
  Methods for Fitted Interactions
}
\description{
  These are methods specifically for the class \code{"fii"} of
  fitted interpoint interactions.
}
\usage{
\method{print}{fii}(x, \dots)

\method{coef}{fii}(object, \dots)

\method{plot}{fii}(x, \dots)

\method{summary}{fii}(object,\dots)

\method{print}{summary.fii}(x, ...)

\method{coef}{summary.fii}(object, ...)
}
\arguments{
  \item{x,object}{
    An object of class \code{"fii"} representing a fitted
    interpoint interaction.
  }
  \item{\dots}{
    Arguments passed to other methods.
  }
}
\details{
  These are methods for the class \code{"fii"}.
  An object of class \code{"fii"} represents a fitted interpoint
  interaction. It is usually obtained by
  using the command \code{\link{fitin}} to extract the fitted
  interaction part of a fitted point process model. 
  See \code{\link{fitin}} for further explanation of this class.
  
  The commands listed here are methods for the generic functions
  \code{\link{print}},
  \code{\link{summary}},
  \code{\link{plot}} 
  and 
  \code{\link{coef}}
  for objects of the class \code{"fii"}.

  Following the usual convention, \code{summary.fii} returns an object of class
  \code{summary.fii}, for which there is a print method.
  The effect is that, when the user types \code{summary(x)},
  the summary is printed, but when the user types \code{y <- summary(x)},
  the summary information is saved.

  The method \code{coef.fii} extracts the canonical coefficients of
  the fitted interaction, and returns them as a numeric vector.
  The method \code{coef.summary.fii} transforms these values into
  quantities that are more easily interpretable, in a format that
  depends on the particular model.
  
  There are also methods for the generic commands \code{\link{reach}} and
  \code{\link{as.interact}}, described elsewhere.
}
\value{
  The \code{print} and \code{plot} methods return \code{NULL}.

  The \code{summary} method returns an object of class
  \code{summary.fii}.

  \code{coef.fii} returns a numeric vector.
  \code{coef.summary.fii} returns data whose structure
  depends on the model.
}
\author{
  \adrian
  
  
}
\seealso{
  \code{\link{fitin}},
  \code{\link{reach.fii}},
  \code{\link{as.interact.fii}}
}
\examples{
  mod <- ppm(cells, ~1, Strauss(0.1))
  f <- fitin(mod)
  f
  summary(f)
  plot(f)
  coef(f)
  coef(summary(f))
}
\keyword{spatial}
\keyword{methods}
