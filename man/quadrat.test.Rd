\name{quadrat.test}
\alias{quadrat.test}
\alias{quadrat.test.ppp}
\alias{quadrat.test.ppm}
\alias{quadrat.test.quadratcount}
\title{Dispersion Test for Spatial Point Pattern Based on
       Quadrat Counts}
\description{
  Performs a test of Complete Spatial Randomness
  for a given point pattern, based on quadrat counts.
  Alternatively performs a goodness-of-fit test of a fitted
  inhomogeneous Poisson model.
  By default performs chi-squared tests; can also perform
  Monte Carlo based tests.
}
\usage{
quadrat.test(X, ...)

\method{quadrat.test}{ppp}(X, nx=5, ny=nx,
                          alternative=c("two.sided", "regular", "clustered"),
                           method=c("Chisq", "MonteCarlo"),
                           conditional=TRUE, CR=1,
                           lambda=NULL, 
                           ...,
                           xbreaks=NULL, ybreaks=NULL, tess=NULL,
                           nsim=1999)

\method{quadrat.test}{ppm}(X, nx=5, ny=nx, 
                          alternative=c("two.sided", "regular", "clustered"),
                           method=c("Chisq", "MonteCarlo"),
                           conditional=TRUE, CR=1, 
                           ..., 
                           xbreaks=NULL, ybreaks=NULL, tess=NULL,
                           nsim=1999)

\method{quadrat.test}{quadratcount}(X,
                          alternative=c("two.sided", "regular", "clustered"),
                          method=c("Chisq", "MonteCarlo"),
                          conditional=TRUE, CR=1,
                          lambda=NULL, 
                          ...,
                          nsim=1999)
}
\arguments{
  \item{X}{
    A point pattern (object of class \code{"ppp"})
    to be subjected to the goodness-of-fit test.
    Alternatively a fitted point process model (object of class
    \code{"ppm"}) to be tested.
    Alternatively \code{X} can be the result of applying
    \code{\link{quadratcount}} to a point pattern.
  }
  \item{nx,ny}{
    Numbers of quadrats in the \eqn{x} and \eqn{y} directions.
    Incompatible with \code{xbreaks} and \code{ybreaks}.
  }
  \item{alternative}{
    Character string (partially matched) specifying the alternative
    hypothesis.
  }
  \item{method}{
    Character string (partially matched) specifying the test to use:
    either \code{method="Chisq"} for the chi-squared test (the default),
    or \code{method="MonteCarlo"} for a Monte Carlo test.
  }
  \item{conditional}{
    Logical.  Should the Monte Carlo test be conducted
    conditionally upon the observed number of points of the pattern?
    Ignored if \code{method="Chisq"}.
  }
  \item{CR}{
    Optional. Numerical value of the index \eqn{\lambda}{lambda}
    for the Cressie-Read test statistic. 
  }
  \item{lambda}{
    Optional. Pixel image (object of class \code{"im"})
    or function (class \code{"funxy"}) giving the predicted
    intensity of the point process.
  }
  \item{\dots}{Ignored.}
  \item{xbreaks}{
    Optional. Numeric vector giving the \eqn{x} coordinates of the
    boundaries of the quadrats. Incompatible with \code{nx}.
  }
  \item{ybreaks}{
    Optional. Numeric vector giving the \eqn{y} coordinates of the
    boundaries of the quadrats. Incompatible with \code{ny}.
  }
  \item{tess}{
    Tessellation (object of class \code{"tess"} or something acceptable
    to \code{\link{as.tess}}) determining the
    quadrats. Incompatible with \code{nx, ny, xbreaks, ybreaks}.
  }
  \item{nsim}{
    The number of simulated samples to generate when
    \code{method="MonteCarlo"}. 
  }
}
\details{
  These functions perform \eqn{\chi^2}{chi^2} tests or Monte Carlo tests
  of goodness-of-fit for a point process model, based on quadrat counts.

  The function \code{quadrat.test} is generic, with methods for
  point patterns (class \code{"ppp"}), split point patterns
  (class \code{"splitppp"}), point process models
  (class \code{"ppm"}) and quadrat count tables (class \code{"quadratcount"}).

  \itemize{
    \item
    if \code{X} is a point pattern, we test the null hypothesis
    that the data pattern is a realisation of Complete Spatial
    Randomness (the uniform Poisson point process). Marks in the point
    pattern are ignored. (If \code{lambda} is given then the null
    hypothesis is the Poisson process with intensity \code{lambda}.)
    \item 
    if \code{X} is a split point pattern, then for each of the
    component point patterns (taken separately) we test 
    the null hypotheses of Complete Spatial Randomness.
    See \code{\link{quadrat.test.splitppp}} for documentation.
    \item
    If \code{X} is a fitted point process model, then it should be
    a Poisson point process model. The 
    data to which this model was fitted are extracted from the model
    object, and are treated as the data point pattern for the test.
    We test the null hypothesis 
    that the data pattern is a realisation of the (inhomogeneous) Poisson point
    process specified by \code{X}.
  }

  In all cases, the window of observation is divided
  into tiles, and the number of data points in each tile is
  counted, as described in \code{\link{quadratcount}}.
  The quadrats are rectangular by default, or may be regions of arbitrary shape
  specified by the argument \code{tess}.
  The expected number of points in each quadrat is also calculated,
  as determined by CSR (in the first case) or by the fitted model
  (in the second case).
  Then the Pearson \eqn{X^2} statistic 
  \deqn{
    X^2 = sum((observed - expected)^2/expected)
  }
  is computed.
  
  If \code{method="Chisq"} then a \eqn{\chi^2}{chi^2} test of
  goodness-of-fit is performed by comparing the test statistic
  to the \eqn{\chi^2}{chi^2} distribution
  with \eqn{m-k} degrees of freedom, where \code{m} is the number of
  quadrats and \eqn{k} is the number of fitted parameters
  (equal to 1 for \code{quadrat.test.ppp}). The default is to
  compute the \emph{two-sided} \eqn{p}-value, so that the test will
  be declared significant if \eqn{X^2} is either very large or very
  small. One-sided \eqn{p}-values can be obtained by specifying the
  \code{alternative}. An important requirement of the
  \eqn{\chi^2}{chi^2} test is that the expected counts in each quadrat
  be greater than 5.

  If \code{method="MonteCarlo"} then a Monte Carlo test is performed,
  obviating the need for all expected counts to be at least 5.  In the
  Monte Carlo test, \code{nsim} random point patterns are generated
  from the null hypothesis (either CSR or the fitted point process
  model). The Pearson \eqn{X^2} statistic is computed as above.
  The \eqn{p}-value is determined by comparing the \eqn{X^2}
  statistic for the observed point pattern, with the values obtained
  from the simulations. Again the default is to
  compute the \emph{two-sided} \eqn{p}-value.

  If \code{conditional} is \code{TRUE} then the simulated samples are
  generated from the multinomial distribution with the number of \dQuote{trials}
  equal to the number of observed points and the vector of probabilities
  equal to the expected counts divided by the sum of the expected counts.
  Otherwise the simulated samples are independent Poisson counts, with
  means equal to the expected counts.

  If the argument \code{CR} is given, then instead of the
  Pearson \eqn{X^2} statistic, the Cressie-Read (1984) power divergence
  test statistic
  \deqn{
    2nI = \frac{2}{\lambda(\lambda+1)}
           \sum_i \left[ \left( \frac{X_i}{E_i} \right)^\lambda - 1 \right]
  }{
    2nI = (2/(lambda * (lambda+1))) * sum((X[i]/E[i])^lambda - 1)
  }
  is computed, where \eqn{X_i}{X[i]} is the \eqn{i}th observed count
  and \eqn{E_i}{E[i]} is the corresponding expected count,
  and the exponent \eqn{\lambda}{lambda} is equal to \code{CR}.
  The value \code{CR=1} gives the Pearson \eqn{X^2} statistic;
  \code{CR=0} gives the likelihood ratio test statistic \eqn{G^2};
  \code{CR=-1/2} gives the Freeman-Tukey statistic \eqn{T^2};
  \code{CR=-1} gives the modified likelihood ratio test statistic \eqn{GM^2};
  and \code{CR=-2} gives Neyman's modified statistic \eqn{NM^2}.
  In all cases the asymptotic distribution of this test statistic is
  the same \eqn{\chi^2}{chi^2} distribution as above.

  The return value is an object of class \code{"htest"}.
  Printing the object gives comprehensible output
  about the outcome of the test.

  The return value also belongs to
  the special class \code{"quadrat.test"}. Plotting the object
  will display the quadrats, annotated by their observed and expected
  counts and the Pearson residuals. See the examples.
}
\seealso{
  \code{\link{quadrat.test.splitppp}},
  \code{\link{quadratcount}},
  \code{\link{quadrats}},
  \code{\link{quadratresample}},
  \code{\link{chisq.test}},
  \code{\link{cdf.test}}.

  To test a Poisson point process model against a specific alternative,
  use \code{\link{anova.ppm}}.
}
\value{
  An object of class \code{"htest"}. See \code{\link{chisq.test}}
  for explanation. 

  The return value is also an object of the special class
  \code{"quadrattest"}, and there is a plot method for this class.
  See the examples.
}
\references{
  Cressie, N. and Read, T.R.C. (1984)
  Multinomial goodness-of-fit tests.
  \emph{Journal of the Royal Statistical Society, Series B}
  \bold{46}, 440--464. 
}
\examples{
  data(simdat)
  quadrat.test(simdat)
  quadrat.test(simdat, 4, 3)

  quadrat.test(simdat, alternative="regular")
  quadrat.test(simdat, alternative="clustered")

  # Using Monte Carlo p-values
  quadrat.test(swedishpines) # Get warning, small expected values.
  \dontrun{
    quadrat.test(swedishpines, method="M", nsim=4999)
    quadrat.test(swedishpines, method="M", nsim=4999, conditional=FALSE)
  }
  \testonly{
    quadrat.test(swedishpines, method="M", nsim=19)
    quadrat.test(swedishpines, method="M", nsim=19, conditional=FALSE)
  }

  # quadrat counts
  qS <- quadratcount(simdat, 4, 3)
  quadrat.test(qS)

  # fitted model: inhomogeneous Poisson
  fitx <- ppm(simdat, ~x, Poisson())
  quadrat.test(fitx)

  te <- quadrat.test(simdat, 4)
  residuals(te)  # Pearson residuals

  plot(te)

  plot(simdat, pch="+", cols="green", lwd=2)
  plot(te, add=TRUE, col="red", cex=1.4, lty=2, lwd=3)

  sublab <- eval(substitute(expression(p[chi^2]==z),
                       list(z=signif(te$p.value,3))))
  title(sub=sublab, cex.sub=3)

  # quadrats of irregular shape
  B <- dirichlet(runifpoint(6, simdat$window))
  qB <- quadrat.test(simdat, tess=B)
  plot(simdat, main="quadrat.test(simdat, tess=B)", pch="+")
  plot(qB, add=TRUE, col="red", lwd=2, cex=1.2)

}
\author{Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  and Rolf Turner \email{r.turner@auckland.ac.nz}
  
}
\keyword{spatial}
\keyword{htest}
