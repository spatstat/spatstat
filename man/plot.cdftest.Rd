\name{plot.cdftest}
\alias{plot.cdftest}
\title{Plot a Spatial Distribution Test}
\description{
  Plot the result of a spatial distribution test
  computed by \code{cdf.test}.
}
\usage{
\method{plot}{cdftest}(x, ...,
                   style=c("cdf", "PP", "QQ"),
                   lwd=par("lwd"), col=par("col"), lty=par("lty"),
                   lwd0=lwd, col0=2, lty0=2,
                   do.legend)
}
\arguments{
  \item{x}{
    Object to be plotted. An object of class \code{"cdftest"}
    produced by a method for \code{\link{cdf.test}}.
  }
  \item{\dots}{
    extra arguments that will be passed to the plotting function
    \code{\link{plot.default}}.
  }
  \item{style}{
    Style of plot. See Details.
  }
  \item{col,lwd,lty}{
    The width, colour and type of lines used to plot the
    empirical curve (the empirical distribution, or PP plot or QQ plot).
  }
  \item{col0,lwd0,lty0}{
    The width, colour and type of lines used to plot the
    reference curve (the predicted distribution, or the diagonal).
  }
  \item{do.legend}{
    Logical value indicating whether to add an
    explanatory legend. Applies only when \code{style="cdf"}.
  }
}
\value{
  \code{NULL}.
}
\details{
  This is the \code{plot} method for the class \code{"cdftest"}.
  An object of this class represents the outcome of
  a spatial distribution test, computed by \code{\link{cdf.test}},
  and based on either the Kolmogorov-Smirnov,
  \ifelse{latex}{\out{Cram\'er}}{Cramer}-von Mises
  or Anderson-Darling test.

  If \code{style="cdf"} (the default), 
  the plot displays the two cumulative distribution functions
  that are compared by the test: namely the empirical cumulative distribution
  function of the covariate at the data points, and the predicted
  cumulative distribution function of the covariate under the model,
  both plotted against the value of the covariate. The
  Kolmogorov-Smirnov test statistic (for example)
  is the maximum vertical separation
  between the two curves.

  If \code{style="PP"} then the P-P plot is drawn. The
  \eqn{x} coordinates of the plot are cumulative
  probabilities for the covariate under the model.
  The \eqn{y} coordinates are cumulative probabilities
  for the covariate at the data points. The diagonal line
  \eqn{y=x} is also drawn for reference. The Kolmogorov-Smirnov
  test statistic is the maximum vertical separation
  between the P-P plot and the diagonal reference line.

  If \code{style="QQ"} then the Q-Q plot is drawn. The
  \eqn{x} coordinates of the plot are quantiles
  of the covariate under the model.
  The \eqn{y} coordinates are quantiles of the 
  covariate at the data points. The diagonal line
  \eqn{y=x} is also drawn for reference. The Kolmogorov-Smirnov
  test statistic cannot be read off the Q-Q plot.
}
\seealso{
  \code{\link{cdf.test}}
}
\examples{
   op <- options(useFancyQuotes=FALSE)

   # synthetic data: nonuniform Poisson process
   X <- rpoispp(function(x,y) { 100 * exp(x) }, win=square(1))

   # fit uniform Poisson process
   fit0 <- ppm(X, ~1)

   # test covariate = x coordinate
   xcoord <- function(x,y) { x }

   # test wrong model
   k <- cdf.test(fit0, xcoord)

   # plot result of test
   plot(k, lwd0=3)

   plot(k, style="PP")

   plot(k, style="QQ")

   options(op)
}

\author{\adrian
  
  
  and \rolf
  
}
\keyword{spatial}
\keyword{hplot}
