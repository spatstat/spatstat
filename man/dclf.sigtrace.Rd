\name{dclf.sigtrace}
\alias{dclf.sigtrace}
\alias{mad.sigtrace}
\alias{mctest.sigtrace}
\title{
  Significance Trace of Cressie-Loosmore-Ford or Maximum Absolute
  Deviation Test
}
\description{
  Generates a Significance Trace of the
  Diggle(1986)/ Cressie (1991)/ Loosmore and Ford (2006) test or the
  Maximum Absolute Deviation test for a spatial point pattern.
}
\usage{
dclf.sigtrace(X, \dots)
mad.sigtrace(X, \dots)
mctest.sigtrace(X, fun=Lest, \dots,
                exponent=1, interpolate=FALSE, alpha=0.05,
                confint=TRUE, rmin=0)
}
\arguments{
  \item{X}{
    Either a point pattern (object of class \code{"ppp"}, \code{"lpp"}
    or other class), a fitted point process model (object of class \code{"ppm"},
    \code{"kppm"} or other class) or an envelope object (class
    \code{"envelope"}). 
  }
  \item{\dots}{
    Arguments passed to \code{\link{envelope}}
    or \code{\link{mctest.progress}}.
    Useful arguments include \code{fun} to determine the summary
    function, \code{nsim} to specify the number of Monte Carlo
    simulations, \code{alternative} to specify a one-sided test,
    and \code{verbose=FALSE} to turn off the messages.
  }
  \item{fun}{
    Function that computes the desired summary statistic
    for a point pattern.
  }
  \item{exponent}{
    Positive number. The exponent of the \eqn{L^p} distance.
    See Details.
  }
  \item{interpolate}{
    Logical value specifying whether to calculate the \eqn{p}-value
    by interpolation.
    If \code{interpolate=FALSE} (the default), a standard Monte Carlo test
    is performed, yielding a \eqn{p}-value of the form \eqn{(k+1)/(n+1)}
    where \eqn{n} is the number of simulations and \eqn{k} is the number
    of simulated values which are more extreme than the observed value.
    If \code{interpolate=TRUE}, the \eqn{p}-value is calculated by
    applying kernel density estimation to the simulated values, and
    computing the tail probability for this estimated distribution.
  }
  \item{alpha}{
    Significance level to be plotted (this has no effect on the calculation
    but is simply plotted as a reference value).
  }
  \item{confint}{
    Logical value indicating whether to compute a confidence interval
    for the \sQuote{true} \eqn{p}-value.
  }
  \item{rmin}{
    Optional. Left endpoint for the interval of \eqn{r} values
    on which the test statistic is calculated.
  }
}
\details{
  The Diggle (1986)/ Cressie (1991)/Loosmore and Ford (2006) test and the 
  Maximum Absolute Deviation test for a spatial point pattern
  are described in \code{\link{dclf.test}}.
  These tests depend on the choice of an interval of
  distance values (the argument \code{rinterval}).
  A \emph{significance trace} (Bowman and Azzalini, 1997;
  Baddeley et al, 2014, 2015)
  of the test is a plot of the \eqn{p}-value
  obtained from the test against the length of
  the interval \code{rinterval}.
  
  The command \code{dclf.sigtrace} performs 
  \code{\link{dclf.test}} on \code{X} using all possible intervals
  of the form \eqn{[0,R]}, and returns the resulting \eqn{p}-values
  as a function of \eqn{R}.

  Similarly \code{mad.sigtrace} performs
  \code{\link{mad.test}} using all possible intervals
  and returns the \eqn{p}-values.

  More generally, \code{mctest.sigtrace} performs a test based on the
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
  
  The result of each command
  is an object of class \code{"fv"} that can be plotted to
  obtain the significance trace. The plot shows the Monte Carlo
  \eqn{p}-value (solid black line), 
  the critical value \code{0.05} (dashed red line),
  and a pointwise 95\% confidence band (grey shading)
  for the \sQuote{true} (Neyman-Pearson) \eqn{p}-value.
  The confidence band is based on the Agresti-Coull (1998)
  confidence interval for a binomial proportion (when
  \code{interpolate=FALSE}) or the delta method
  and normal approximation (when \code{interpolate=TRUE}).

  If \code{X} is an envelope object and \code{fun=NULL} then
  the code will re-use the simulated functions stored in \code{X}.
}
\value{
  An object of class \code{"fv"} that can be plotted to
  obtain the significance trace. 
}
\references{
  Agresti, A. and Coull, B.A. (1998)
  Approximate is better than \dQuote{Exact} for interval
   estimation of binomial proportions.
  \emph{American Statistician} \bold{52}, 119--126.

  Baddeley, A., Diggle, P., Hardegen, A., Lawrence, T.,
  Milne, R. and Nair, G. (2014)
  On tests of spatial pattern based on simulation envelopes.
  \emph{Ecological Monographs} \bold{84}(3) 477--489.

  Baddeley, A., Hardegen, A., Lawrence, L., 
  Milne, R.K., Nair, G.M. and Rakshit, S. (2015)
  Pushing the envelope: extensions of graphical
  Monte Carlo tests. Submitted for publication.

  Bowman, A.W. and Azzalini, A. (1997) 
  \emph{Applied smoothing techniques for data analysis: 
    the kernel approach with S-Plus illustrations}.
  Oxford University Press, Oxford.
}
\author{
  \adrian, Andrew Hardegen, Tom Lawrence,
  Robin Milne, Gopalan Nair and Suman Rakshit.
  Implemented by
  \adrian
  
  
  \rolf
  
  and \ege
  
}
\seealso{
  \code{\link{dclf.test}} for the tests;
  \code{\link{dclf.progress}} for progress plots.
  
  See \code{\link{plot.fv}} for information on plotting
  objects of class \code{"fv"}.

  See also \code{\link{dg.sigtrace}}.
}
\examples{
  plot(dclf.sigtrace(cells, Lest, nsim=19))
}
\keyword{spatial}
\keyword{htest}
