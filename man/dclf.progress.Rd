\name{dclf.progress}
\alias{dclf.progress}
\alias{mad.progress}
\alias{mctest.progress}
\title{
  Progress Plot of Test of Spatial Pattern
}
\description{
  Generates a progress plot (envelope representation) of the
  Diggle-Cressie-Loosmore-Ford test or the
  Maximum Absolute Deviation test for a spatial point pattern.
}
\usage{
dclf.progress(X, \dots)
mad.progress(X, \dots)
mctest.progress(X, fun = Lest, \dots,
                exponent = 1, nrank = 1,
                interpolate = FALSE, alpha, rmin=0)
}
\arguments{
  \item{X}{
    Either a point pattern (object of class \code{"ppp"}, \code{"lpp"}
    or other class), a fitted point process model (object of class \code{"ppm"},
    \code{"kppm"} or other class) or an envelope object (class
    \code{"envelope"}). 
  }
  \item{\dots}{
    Arguments passed to \code{mctest.progress} or to \code{\link{envelope}}.
    Useful arguments include \code{fun} to determine the summary
    function, \code{nsim} to specify the number of Monte Carlo
    simulations, \code{alternative} to specify one-sided or two-sided
    envelopes, and \code{verbose=FALSE} to turn off the messages.
  }
  \item{fun}{
    Function that computes the desired summary statistic
    for a point pattern. 
  }
  \item{exponent}{
    Positive number. The exponent of the \eqn{L^p} distance.
    See Details.
  }
  \item{nrank}{
    Integer. The rank of the critical value of the Monte Carlo test,
    amongst the \code{nsim} simulated values.
    A rank of 1 means that the minimum and maximum
    simulated values will become the critical values for the test.
  }
  \item{interpolate}{
    Logical value indicating how to compute the critical value.
    If \code{interpolate=FALSE} (the default), a standard Monte Carlo test
    is performed, and the critical value is the largest
    simulated value of the test statistic (if \code{nrank=1})
    or the \code{nrank}-th largest (if \code{nrank} is another number).
    If \code{interpolate=TRUE}, kernel density estimation
    is applied to the simulated values, and the critical value is
    the upper \code{alpha} quantile of this estimated distribution.
  }
  \item{alpha}{
    Optional. The significance level of the test.
    Equivalent to \code{nrank/(nsim+1)} where \code{nsim} is the
    number of simulations.
  }
  \item{rmin}{
    Optional. Left endpoint for the interval of \eqn{r} values
    on which the test statistic is calculated.
  }
}
\details{
  The Diggle-Cressie-Loosmore-Ford test and the 
  Maximum Absolute Deviation test for a spatial point pattern
  are described in \code{\link{dclf.test}}.
  These tests depend on the choice of an interval of
  distance values (the argument \code{rinterval}).
  A \emph{progress plot} or \emph{envelope representation}
  of the test (Baddeley et al, 2014) is a plot of the
  test statistic (and the corresponding critical value) against the length of
  the interval \code{rinterval}.
  
  The command \code{dclf.progress} performs 
  \code{\link{dclf.test}} on \code{X} using all possible intervals
  of the form \eqn{[0,R]}, and returns the resulting values of the test
  statistic, and the corresponding critical values of the test,
  as a function of \eqn{R}. 

  Similarly \code{mad.progress} performs
  \code{\link{mad.test}} using all possible intervals
  and returns the test statistic and critical value.

  More generally, \code{mctest.progress} performs a test based on the
  \eqn{L^p} discrepancy between the curves. The deviation between two
  curves is measured by the \eqn{p}th root of the integral of
  the \eqn{p}th power of the absolute value of the difference
  between the two curves. The exponent \eqn{p} is
  given by the argument \code{exponent}. The case \code{exponent=2}
  is the Cressie-Loosmore-Ford test, while \code{exponent=Inf} is the
  MAD test.

  If the argument \code{rmin} is given, it specifies the left endpoint
  of the interval defining the test statistic: the tests are
  performed using intervals \eqn{[r_{\mbox{\scriptsize min}},R]}{[rmin,R]}
  where \eqn{R \ge r_{\mbox{\scriptsize min}}}{R >= rmin}.
  
  The result of each command is an object of class \code{"fv"}
  that can be plotted to obtain the progress plot. The display shows
  the test statistic (solid black line) and the Monte Carlo
  acceptance region (grey shading).

  The significance level for the Monte Carlo test is
  \code{nrank/(nsim+1)}. Note that \code{nsim} defaults to 99,
  so if the values of \code{nrank} and \code{nsim} are not given,
  the default is a test with significance level 0.01.

  If \code{X} is an envelope object, then some of the data stored
  in \code{X} may be re-used:
  \itemize{
    \item
    If \code{X} is an envelope object containing simulated functions,
    and \code{fun=NULL}, then
    the code will re-use the simulated functions stored in \code{X}.
    \item
    If \code{X} is an envelope object containing
    simulated point patterns, 
    then \code{fun} will be applied to the stored point patterns
    to obtain the simulated functions.
    If \code{fun} is not specified, it defaults to \code{\link{Lest}}.
    \item
    Otherwise, new simulations will be performed,
    and \code{fun} defaults to  \code{\link{Lest}}.
  }
}
\value{
  An object of class \code{"fv"} that can be plotted to
  obtain the progress plot. 
}
\references{
  Baddeley, A., Diggle, P., Hardegen, A., Lawrence, T.,
  Milne, R. and Nair, G. (2014)
  On tests of spatial pattern based on simulation envelopes.
  \emph{Ecological Monographs} \bold{84} (3) 477--489.
}
\author{
  \adrian
  
  ,
  Andrew Hardegen, Tom Lawrence, Gopal Nair and
  Robin Milne.
}
\seealso{
  \code{\link{dclf.test}} and
  \code{\link{mad.test}} for the tests.
  
  See \code{\link{plot.fv}} for information on plotting
  objects of class \code{"fv"}.
}
\examples{
  plot(dclf.progress(cells, nsim=19))
}
\keyword{spatial}
\keyword{htest}

