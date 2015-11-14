\name{cdf.test} 
\alias{cdf.test}
\alias{cdf.test.ppm}
\alias{cdf.test.lppm}
\alias{cdf.test.lpp}
\alias{cdf.test.ppp}
\alias{cdf.test.slrm}
\title{Spatial Distribution Test for Point Pattern or Point Process Model} 
\description{
  Performs a test of goodness-of-fit of a point process model.
  The observed and predicted distributions
  of the values of a spatial covariate are compared using either the
  Kolmogorov-Smirnov test,
  \ifelse{latex}{\out{Cram\'er}}{Cramer}-von Mises test
  or Anderson-Darling test.
  For non-Poisson models, a Monte Carlo test is used.
}
\usage{
cdf.test(...)

\method{cdf.test}{ppp}(X, covariate, test=c("ks", "cvm", "ad"), \dots, jitter=TRUE)

\method{cdf.test}{ppm}(model, covariate,  test=c("ks", "cvm", "ad"), \dots,
          jitter=TRUE, nsim=99, verbose=TRUE)

\method{cdf.test}{lpp}(X, covariate,  test=c("ks", "cvm", "ad"), \dots, jitter=TRUE)

\method{cdf.test}{lppm}(model, covariate,  test=c("ks", "cvm", "ad"),
\dots,
          jitter=TRUE, nsim=99, verbose=TRUE)

\method{cdf.test}{slrm}(model, covariate,  test=c("ks", "cvm", "ad"), \dots, modelname=NULL, covname=NULL)
}
\arguments{
  \item{X}{
    A point pattern (object of class \code{"ppp"} or \code{"lpp"}).
  }
  \item{model}{
    A fitted point process model (object of class \code{"ppm"} or \code{"lppm"})
    or fitted spatial logistic regression (object of class \code{"slrm"}).
  }
  \item{covariate}{
    The spatial covariate on which the test will be based.
    A function, a pixel image (object of class \code{"im"}),
    a list of pixel images, or one of the characters
    \code{"x"} or \code{"y"} indicating the Cartesian coordinates.
  }
  \item{test}{
    Character string identifying the test to be performed:
    \code{"ks"} for Kolmogorov-Smirnov test,
    \code{"cvm"} for \ifelse{latex}{\out{Cram\'er}}{Cramer}-von Mises test
    or \code{"ad"} for Anderson-Darling test.
  }
  \item{\dots}{
    Arguments passed to \code{\link[stats]{ks.test}}
    (from the \pkg{stats} package) or 
    \code{\link[goftest]{cvm.test}} or
    \code{\link[goftest]{ad.test}} (from the \pkg{goftest} package)
    to control the test.
  }
  \item{jitter}{
    Logical flag. If \code{jitter=TRUE}, values of the covariate
    will be slightly perturbed at random, to avoid tied values in the test.
  }
  \item{modelname,covname}{
    Character strings giving alternative names for \code{model}
    and \code{covariate} to be used in labelling plot axes.
  }
  \item{nsim}{
    Number of simulated realisations from the \code{model} to be used
    for the Monte Carlo test, when \code{model} is not a Poisson process.
  }
  \item{verbose}{
    Logical value indicating whether to print progress reports when
    performing a Monte Carlo test.
  }
}
\details{
  These functions perform a goodness-of-fit test of a Poisson or Gibbs point
  process model fitted to point pattern data. The observed distribution
  of the values of a spatial covariate at the data points,
  and the predicted distribution of the same values under the model,
  are compared using the Kolmogorov-Smirnov test,
  the \ifelse{latex}{\out{Cram\'er}}{Cramer}-von Mises test
  or the Anderson-Darling test. For Gibbs models, a Monte Carlo test is
  performed using these test statistics.

  The function \code{cdf.test} is generic, with methods for
  point patterns (\code{"ppp"} or \code{"lpp"}),
  point process models (\code{"ppm"} or \code{"lppm"})
  and spatial logistic regression models (\code{"slrm"}).
  \itemize{
    \item 
    If \code{X} is a point pattern dataset (object of class
    \code{"ppp"}), then \code{cdf.test(X, \dots)}
    performs a goodness-of-fit test of the
    uniform Poisson point process (Complete Spatial Randomness, CSR)
    for this dataset.
    For a multitype point pattern, the uniform intensity
    is assumed to depend on the type of point (sometimes called
    Complete Spatial Randomness and Independence, CSRI).
    \item
    If \code{model} is a fitted point process model
    (object of class \code{"ppm"} or \code{"lppm"})
    then \code{cdf.test(model, \dots)} performs
    a test of goodness-of-fit for this fitted model. 
    \item
    If \code{model} is a fitted spatial logistic regression
    (object of class \code{"slrm"}) then \code{cdf.test(model, \dots)} performs
    a test of goodness-of-fit for this fitted model. 
  }
  The test is performed by comparing the observed distribution
  of the values of a spatial covariate at the data points,
  and the predicted distribution of the same covariate under the model,
  using a classical goodness-of-fit test.
  Thus, you must nominate
  a spatial covariate for this test.
  
  If \code{X} is a point pattern that does not have marks,
  the argument \code{covariate} should be either a \code{function(x,y)}
  or a pixel image (object of class \code{"im"} containing the values
  of a spatial function, or one of the characters \code{"x"} or
  \code{"y"} indicating the Cartesian coordinates.
  If \code{covariate} is an image, it should have numeric values,
  and its domain should cover the observation window of the
  \code{model}. If \code{covariate} is a function, it should expect
  two arguments \code{x} and \code{y} which are vectors of coordinates,
  and it should return a numeric vector of the same length
  as \code{x} and \code{y}.
  
  If \code{X} is a multitype point pattern, the argument \code{covariate}
  can be either a \code{function(x,y,marks)},
  or a pixel image, or a list of pixel images corresponding to
  each possible mark value, or one of the characters \code{"x"} or
  \code{"y"} indicating the Cartesian coordinates. 
  
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
  these numbers are i.i.d. uniform random numbers. The
  A goodness-of-fit test of the uniform distribution is applied
  to these numbers using \code{stats::\link[stats]{ks.test}},
  \code{goftest::\link[goftest]{cvm.test}} or
  \code{goftest::\link[goftest]{ad.test}}.

  This test was apparently first described (in the context of
  spatial data, and using Kolmogorov-Smirnov) by Berman (1986).
  See also Baddeley et al (2005).

  If \code{model} is not a Poisson process, then
  a Monte Carlo test is performed, by generating \code{nsim}
  point patterns which are simulated realisations of the \code{model},
  re-fitting the model to each simulated point pattern, 
  and calculating the test statistic for each fitted model.
  The Monte Carlo \eqn{p} value is determined by comparing
  the simulated values of the test statistic 
  with the value for the original data.
  
  The return value is an object of class \code{"htest"} containing the
  results of the hypothesis test. The print method for this class
  gives an informative summary of the test outcome.

  The return value also belongs to the class \code{"cdftest"}
  for which there is a plot method \code{\link{plot.cdftest}}.
  The plot method displays the empirical cumulative distribution
  function of the covariate at the data points, and the predicted
  cumulative distribution function of the covariate under the model,
  plotted against the value of the covariate.

  The argument \code{jitter} controls whether covariate values are
  randomly perturbed, in order to avoid ties.
  If the original data contains any ties in the covariate (i.e. points
  with equal values of the covariate), and if \code{jitter=FALSE}, then 
  the Kolmogorov-Smirnov test implemented in \code{\link[stats]{ks.test}}
  will issue a warning that it cannot calculate the exact \eqn{p}-value.
  To avoid this, if \code{jitter=TRUE} each value of the covariate will
  be perturbed by adding a small random value. The perturbations are
  normally distributed with standard deviation equal to one hundredth of
  the range of values of the covariate. This prevents ties, 
  and the \eqn{p}-value is still correct. There is
  a very slight loss of power.
}
\value{
  An object of class \code{"htest"} containing the results of the
  test. See \code{\link[stats]{ks.test}} for details. The return value can be
  printed to give an informative summary of the test.

  The value also belongs to the class \code{"cdftest"} for which there is
  a plot method.
}
\section{Warning}{
  The outcome of the test involves a small amount of random variability,
  because (by default) the coordinates are randomly perturbed to
  avoid tied values. Hence, if \code{cdf.test} is executed twice, the
  \eqn{p}-values will not be exactly the same. To avoid this behaviour,
  set \code{jitter=FALSE}.
}
\author{\adrian
  
  
  and \rolf
  
}
\seealso{
  \code{\link{plot.cdftest}},
  \code{\link{quadrat.test}},
  \code{\link{berman.test}},
  \code{\link[stats]{ks.test}},
  \code{\link[goftest]{cvm.test}},
  \code{\link[goftest]{ad.test}},
  \code{\link{ppm}}
}
\references{
  Baddeley, A., Turner, R.,
  \ifelse{latex}{\out{M\o ller}}{Moller}, J. and Hazelton, M. (2005)
  Residual analysis for spatial point processes.
  \emph{Journal of the Royal Statistical Society, Series B}
  \bold{67}, 617--666.

  Berman, M. (1986)
  Testing for spatial association between a point process
  and another stochastic process.
  \emph{Applied Statistics} \bold{35}, 54--62.
}
\examples{
   op <- options(useFancyQuotes=FALSE)

   # test of CSR using x coordinate
   cdf.test(nztrees, "x")
   cdf.test(nztrees, "x", "cvm")
   cdf.test(nztrees, "x", "ad")

   # test of CSR using a function of x and y
   fun <- function(x,y){2* x + y}
   cdf.test(nztrees, fun)

   # test of CSR using an image covariate
   funimage <- as.im(fun, W=Window(nztrees))
   cdf.test(nztrees, funimage)

   # fit inhomogeneous Poisson model and test
   model <- ppm(nztrees ~x)
   cdf.test(model, "x")

   if(interactive()) {
     # synthetic data: nonuniform Poisson process
     X <- rpoispp(function(x,y) { 100 * exp(x) }, win=square(1))

     # fit uniform Poisson process
     fit0 <- ppm(X ~1)
     # fit correct nonuniform Poisson process
     fit1 <- ppm(X ~x)

     # test wrong model
     cdf.test(fit0, "x")
     # test right model
     cdf.test(fit1, "x")
   }

   # multitype point pattern
   cdf.test(amacrine, "x")
   yimage <- as.im(function(x,y){y}, W=Window(amacrine))
   cdf.test(ppm(amacrine ~marks+y), yimage)

   options(op)
}
\keyword{htest}
\keyword{spatial}

