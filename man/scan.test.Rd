\name{scan.test}
\alias{scan.test}
\title{
  Spatial Scan Test
}
\description{
  Performs the Spatial Scan Test for clustering
  in a spatial point pattern, or for clustering of one type of point
  in a bivariate spatial point pattern.
}
\usage{
scan.test(X, r, ...,
          method = c("poisson", "binomial"),
          nsim = 19,
          baseline = NULL,
          case = 2, 
          alternative = c("greater", "less", "two.sided"),
          verbose = TRUE)
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
  \item{nsim}{
    Number of simulations for computing Monte Carlo p-value.
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
  \item{verbose}{
    Logical. Whether to print progress reports.
  }
}
\details{
  The spatial scan test (Kulldorf, 1997) is applied
  to the point pattern \code{X}.

  In a nutshell,
  \itemize{
    \item
    If \code{method="poisson"} then 
    a significant result would mean that there is a circle of radius
    \code{r}, located somewhere in the spatial domain of the data,
    which contains a significantly higher than
    expected number of points of \code{X}. That is, the
    pattern \code{X} exhibits spatial clustering.
    \item 
    If \code{method="binomial"} then \code{X} must be a bivariate (two-type)
    point pattern. By default, the first type of point is interpreted as
    a control (non-event) and the second type of point as a case (event).
    A significant result would mean that there is a
    circle of radius \code{r} which contains a significantly higher than
    expected number of cases. That is, the cases are clustered together,
    conditional on the locations of all points.
  }

  Following is a more detailed explanation.
  \itemize{
    \item 
    If \code{method="poisson"} then the scan test based on Poisson
    likelihood is performed (Kulldorf, 1997).
    The dataset \code{X} is treated as an unmarked point pattern.
    By default (if \code{baseline} is not specified) 
    the null hypothesis is complete spatial randomness CSR
    (i.e. a uniform Poisson process).
    The alternative hypothesis is a Poisson process with
    one intensity \eqn{\beta_1}{beta1} inside some circle of radius
    \code{r} and another intensity \eqn{\beta_0}{beta0} outside the
    circle.
    If \code{baseline} is given, then it should be a pixel image
    or a \code{function(x,y)}. The null hypothesis is
    an inhomogeneous Poisson process with intensity proportional
    to \code{baseline}. The alternative hypothesis is an inhomogeneous
    Poisson process with intensity
    \code{beta1 * baseline} inside some circle of radius \code{r},
    and \code{beta0 * baseline} outside the circle.
    \item
    If \code{method="binomial"} then the scan test based on
    binomial likelihood is performed (Kulldorf, 1997).
    The dataset \code{X} must be a bivariate point pattern,
    i.e. a multitype point pattern with two types.
    The null hypothesis is that all permutations of the type labels are
    equally likely.
    The alternative hypothesis is that some circle of radius
    \code{r} has a higher proportion of points of the second type,
    than expected under the null hypothesis.
  }

  The result of \code{scan.test} is a hypothesis test
  (object of class \code{"htest"}) which can be plotted to
  report the results. The component \code{p.value} contains the
  \eqn{p}-value.

  The result of \code{scan.test} can also be plotted (using the plot
  method for the class \code{"scan.test"}). The plot is
  a pixel image of the Likelihood Ratio Test Statistic
  (2 times the log likelihood ratio) as a function
  of the location of the centre of the circle.
  This pixel image can be extracted from the object
  using \code{\link{as.im.scan.test}}.
  The Likelihood Ratio Test Statistic is computed by
  \code{\link{scanLRTS}}.
}
\value{
  An object of class \code{"htest"} (hypothesis test)
  which also belongs to the class \code{"scan.test"}.
  Printing this object gives the result of the test.
  Plotting this object displays the Likelihood Ratio Test Statistic
  as a function of the location of the centre of the circle.
}
\references{
  Kulldorff, M. (1997)
  A spatial scan statistic.
  \emph{Communications in Statistics --- Theory and Methods}
  \bold{26}, 1481--1496.
}
\author{Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  and Rolf Turner \email{r.turner@auckland.ac.nz}
  
}
\seealso{
  \code{\link{plot.scan.test}},
  \code{\link{as.im.scan.test}},
  \code{\link{relrisk}},
  \code{\link{scanLRTS}}
}
\examples{
   nsim <- if(interactive()) 19 else 2
   rr <- if(interactive()) seq(0.5, 1, by=0.1) else c(0.5, 1)
   scan.test(redwood, 0.1 * rr, method="poisson", nsim=nsim)
   scan.test(chorley, rr, method="binomial", case="larynx", nsim=nsim)
}
\keyword{htest}
\keyword{spatial}
