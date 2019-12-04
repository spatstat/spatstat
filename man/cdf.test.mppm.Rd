\name{cdf.test.mppm}
\alias{cdf.test.mppm}
\title{Spatial Distribution Test for Multiple Point Process Model} 
\description{
  Performs a spatial distribution test
  of a point process model fitted to multiple spatial point
  patterns.  The test compares the observed
  and predicted distributions of the values of a spatial covariate,
  using either the Kolmogorov-Smirnov,
  \ifelse{latex}{\out{Cram\'er}}{Cramer}-von Mises 
  or Anderson-Darling test of goodness-of-fit. 
}
\usage{
\method{cdf.test}{mppm}(model, covariate, test=c("ks", "cvm", "ad"), ...,
            nsim=19, verbose=TRUE, interpolate=FALSE, fast=TRUE, jitter=TRUE)
}
\arguments{
  \item{model}{
    An object of class \code{"mppm"} representing a point process model
    fitted to multiple spatial point patterns.
  }
  \item{covariate}{
    The spatial covariate on which the test will be based.
    A function, a pixel image, a list of functions, a list of pixel
    images, a hyperframe, a character string containing the name
    of one of the covariates in \code{model}, or one of the strings
    \code{"x"} or \code{"y"}.
  }
  \item{test}{
    Character string identifying the test to be performed:
    \code{"ks"} for Kolmogorov-Smirnov test,
    \code{"cvm"} for \ifelse{latex}{\out{Cram\'er}}{Cramer}-von Mises test
    or \code{"ad"} for Anderson-Darling test.
  }
  \item{\dots}{
    Arguments passed to \code{\link{cdf.test}}
    to control the test.
  }
  \item{nsim}{
    Number of simulated realisations which should be generated,
    if a Monte Carlo test is required.
  }
  \item{verbose}{Logical flag indicating whether to print
    progress reports.
  }
  \item{interpolate}{
    Logical flag indicating whether to interpolate between
    pixel values when code{covariate} is a pixel image.
    See \emph{Details}.
  }
  \item{fast}{
    Logical flag. If \code{TRUE}, values of the covariate
    are only sampled at the original quadrature points used to
    fit the model. If \code{FALSE}, values of the covariate
    are sampled at all pixels, which can be slower by three orders of
    magnitude.
  }
  \item{jitter}{
    Logical flag. If \code{TRUE}, observed values of the covariate
    are perturbed by adding small random values, to avoid
    tied observations.
  }
}
\details{
  This function is a method for the generic function
  \code{\link{cdf.test}} for the class \code{mppm}.

  This function performs a goodness-of-fit test of
  a point process model that has been fitted to multiple point patterns.
  The observed distribution
  of the values of a spatial covariate at the data points,
  and the predicted distribution of the same values under the model,
  are compared using the Kolmogorov-Smirnov,
  \ifelse{latex}{\out{Cram\'er}}{Cramer}-von Mises 
  or Anderson-Darling test of goodness-of-fit.
  These are exact tests if the model is Poisson;
  otherwise, for a Gibbs model, a Monte Carlo p-value is computed by
  generating simulated realisations of the model and applying the
  selected goodness-of-fit test to each simulation.
  
  The argument \code{model} should be a fitted point process model
  fitted to multiple point patterns
  (object of class \code{"mppm"}). 
  
  The argument \code{covariate} contains the values of a spatial
  function. It can be
  \itemize{
    \item a \code{function(x,y)}
    \item a pixel image (object of class \code{"im"}
    \item a list of \code{function(x,y)}, one for each point pattern
    \item a list of pixel images, one for each point pattern
    \item a hyperframe (see \code{\link{hyperframe}})
    of which the first
    column will be taken as containing the covariate
    \item a character string giving the name of one of the covariates
    in \code{model}
    \item one of the character strings \code{"x"} or \code{"y"},
    indicating the spatial coordinates.
  }
  If \code{covariate} is an image, it should have numeric values,
  and its domain should cover the observation window of the
  \code{model}. If \code{covariate} is a function, it should expect
  two arguments \code{x} and \code{y} which are vectors of coordinates,
  and it should return a numeric vector of the same length
  as \code{x} and \code{y}.  

  First the original data point pattern is extracted from \code{model}.
  The values of the \code{covariate} at these data points are
  collected. 

  The predicted distribution of the values of the \code{covariate}
  under the fitted \code{model} is computed as follows.
  The values of the \code{covariate} at all locations in the
  observation window are evaluated,
  weighted according to the point process intensity of the fitted model,
  and compiled into a cumulative distribution function \eqn{F} using
  \code{\link{ewcdf}}.

  The probability integral transformation is then applied:
  the values of the \code{covariate} at the original data points
  are transformed by the predicted cumulative distribution function
  \eqn{F} into numbers between 0 and 1. If the model is correct,
  these numbers are i.i.d. uniform random numbers. 
  A goodness-of-fit test of the uniform distribution is applied
  to these numbers using \code{\link[stats]{ks.test}},
  \code{\link[goftest]{cvm.test}} or \code{\link[goftest]{ad.test}}.

  The argument \code{interpolate} determines 
  how pixel values will be handled when code{covariate} is a pixel image.
  The value of the covariate at a data point is obtained
  by looking up the value of the nearest pixel if
  \code{interpolate=FALSE}, or by linearly interpolating
  between the values of the four nearest pixels 
  if \code{interpolate=TRUE}. Linear interpolation is slower,
  but is sometimes necessary to avoid tied values of the covariate
  arising when the pixel grid is coarse.

  If \code{model} is a Poisson point process, then the 
  Kolmogorov-Smirnov,
  \ifelse{latex}{\out{Cram\'er}}{Cramer}-von Mises 
  and Anderson-Darling tests are theoretically exact.
  This test was apparently first described (in the context of
  spatial data, and for Kolmogorov-Smirnov) by Berman (1986).
  See also Baddeley et al (2005).

  If \code{model} is not a Poisson point process, then the
  Kolmogorov-Smirnov,
  \ifelse{latex}{\out{Cram\'er}}{Cramer}-von Mises 
  and Anderson-Darling tests are biased.
  Instead they are used as the basis of a Monte Carlo test.
  First \code{nsim} simulated realisations of the model will be generated.
  Each simulated realisation consists of a list of simulated point
  patterns, one for each of the original data patterns. This
  can take a very long time. The model is then re-fitted to each
  simulation, and the refitted model is subjected to the goodness-of-fit
  test described above. A Monte Carlo p-value is then computed by
  comparing the p-value of the original test with the
  p-values obtained from the simulations.
}
\value{
  An object of class \code{"cdftest"} and \code{"htest"}
  containing the results of the
  test. See \code{\link{cdf.test}} for details.
}
\seealso{
  \code{\link{cdf.test}},
  \code{\link{quadrat.test}},
  \code{\link{mppm}}
}
\references{
  Baddeley, A., Rubak, E. and Turner, R. (2015)
  \emph{Spatial Point Patterns: Methodology and Applications with R}.
  London: Chapman and Hall/CRC Press. 

  Baddeley, A., Turner, R., \Moller, J. and Hazelton, M. (2005)
  Residual analysis for spatial point processes.
  \emph{Journal of the Royal Statistical Society, Series B}
  \bold{67}, 617--666.

  Berman, M. (1986)
  Testing for spatial association between a point process
  and another stochastic process.
  \emph{Applied Statistics} \bold{35}, 54--62.
}
\author{
  \adrian, Ida-Maria Sintorn and Leanne Bischoff.
  Implemented by 
  \spatstatAuthors.
}
\examples{
   # three i.i.d. realisations of nonuniform Poisson process
   lambda <- as.im(function(x,y) { 300 * exp(x) }, square(1))
   dat <- hyperframe(X=list(rpoispp(lambda), rpoispp(lambda), rpoispp(lambda)))

   # fit uniform Poisson process
   fit0 <- mppm(X~1, dat)
   # fit correct nonuniform Poisson process
   fit1 <- mppm(X~x, dat)

   # test wrong model
   cdf.test(fit0, "x")
   # test right model
   cdf.test(fit1, "x")

   # Gibbs model
   fitGibbs <- update(fit0, interaction=Strauss(0.07))
   ns <- if(interactive()) 19 else 3
   cdf.test(fitGibbs, "x", nsim=ns)
}
\keyword{htest}
\keyword{spatial}

