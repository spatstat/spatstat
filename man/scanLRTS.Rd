\name{scanLRTS}
\alias{scanLRTS}
\title{
  Likelihood Ratio Test Statistic for Scan Test
}
\description{
  Calculate the Likelihood Ratio Test Statistic for the Scan Test,
  at each spatial location.
}
\usage{
scanLRTS(X, r, \dots,
   method = c("poisson", "binomial"),
   baseline = NULL, case = 2,
   alternative = c("greater", "less", "two.sided"),
   saveopt = FALSE,
   Xmask = NULL)
}
\arguments{
  \item{X}{
    A point pattern (object of class \code{"ppp"}).
  }
  \item{r}{
    Radius of circle to use. A single number or a numeric vector.
  }
  \item{\dots}{
    Optional. Arguments passed to \code{\link{as.mask}} to determine the
    spatial resolution of the computations.
  }
  \item{method}{
    Either \code{"poisson"} or \code{"binomial"}
    specifying the type of likelihood.
  }
  \item{baseline}{
    Baseline for the Poisson intensity, if \code{method="poisson"}.
    A pixel image or a function.
  }
  \item{case}{
    Which type of point should be interpreted as a case,
    if \code{method="binomial"}.
    Integer or character string.
  }
  \item{alternative}{
    Alternative hypothesis: \code{"greater"} if the alternative
    postulates that the mean number of points inside the circle
    will be greater than expected under the null.
  }
  \item{saveopt}{
    Logical value indicating to save the optimal value of \code{r}
    at each location.
  }
  \item{Xmask}{
    Internal use only.
  }
}
\details{
  This command computes, for all spatial locations \code{u},
  the Likelihood Ratio Test Statistic \eqn{\Lambda(u)}{Lambda(u)}
  for a test of homogeneity at the location \eqn{u}, as described
  below. The result is a pixel image giving the values of
  \eqn{\Lambda(u)}{Lambda(u)} at each pixel. 

  The \bold{maximum} value of \eqn{\Lambda(u)}{Lambda(u)} over all locations
  \eqn{u} is the \emph{scan statistic}, which is the basis of
  the   \emph{scan test} performed by \code{\link{scan.test}}.

  \itemize{
    \item 
    If \code{method="poisson"} then the test statistic is based on Poisson
    likelihood.
    The dataset \code{X} is treated as an unmarked point pattern.
    By default (if \code{baseline} is not specified) 
    the null hypothesis is complete spatial randomness CSR
    (i.e. a uniform Poisson process).
    At the spatial location \eqn{u},
    the alternative hypothesis is a Poisson process with
    one intensity \eqn{\beta_1}{beta1} inside the circle of radius
    \code{r} centred at \eqn{u},
    and another intensity \eqn{\beta_0}{beta0} outside the
    circle.
    If \code{baseline} is given, then it should be a pixel image
    or a \code{function(x,y)}. The null hypothesis is
    an inhomogeneous Poisson process with intensity proportional
    to \code{baseline}. The alternative hypothesis is an inhomogeneous
    Poisson process with intensity
    \code{beta1 * baseline} inside the circle,
    and \code{beta0 * baseline} outside the circle.
    \item
    If \code{method="binomial"} then the test statistic is based on
    binomial likelihood.
    The dataset \code{X} must be a bivariate point pattern,
    i.e. a multitype point pattern with two types.
    The null hypothesis is that all permutations of the type labels are
    equally likely.
    The alternative hypothesis is that the circle of radius
    \code{r} centred at \eqn{u}
    has a higher proportion of points of the second type,
    than expected under the null hypothesis.
  }

  If \code{r} is a vector of more than one value for the radius,
  then the calculations described above are performed for
  every value of \code{r}. Then the maximum over \code{r} is taken
  for each spatial location \eqn{u}.
  The resulting pixel value of \code{scanLRTS} at a location
  \eqn{u} is the profile maximum of the Likelihood Ratio Test Statistic,
  that is, the maximum of the
  Likelihood Ratio Test Statistic for circles of all radii,
  centred at the same location \eqn{u}.

  If you have already performed a scan test using
  \code{\link{scan.test}}, the  Likelihood Ratio Test Statistic
  can be extracted from the test result using the 
  function \code{\link{as.im.scan.test}}.
}
\section{Warning: window size}{
  Note that the result of \code{scanLRTS} is a pixel image
  on a larger window than the original window of \code{X}.
  The expanded window contains the centre of any circle
  of radius \code{r}
  that has nonempty intersection with the original window.
}
\value{
  A pixel image (object of class \code{"im"}) whose pixel values
  are the values of the (profile) Likelihood Ratio Test Statistic at each
  spatial location. 
}
\references{
  Kulldorff, M. (1997)
  A spatial scan statistic.
  \emph{Communications in Statistics --- Theory and Methods}
  \bold{26}, 1481--1496.
}
\author{\adrian
  
  
  and \rolf
  
}
\seealso{
 \code{\link{scan.test}},
 \code{\link{as.im.scan.test}}
}
\examples{
   plot(scanLRTS(redwood, 0.1, method="poisson"))
   sc <- scanLRTS(chorley, 1, method="binomial", case="larynx") 
   plot(sc)
   scanstatchorley <- max(sc)
}
\keyword{htest}
\keyword{spatial}
