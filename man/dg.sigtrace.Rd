\name{dg.sigtrace}
\alias{dg.sigtrace}
\title{
  Significance Trace of Dao-Genton Test
}
\description{
  Generates a Significance Trace of the
  Dao and Genton (2014) test for a spatial point pattern.
}
\usage{
  dg.sigtrace(X, fun = Lest, \dots,
              exponent = 2, nsim = 19, nsimsub = nsim - 1,
              alternative = c("two.sided", "less", "greater"),
              rmin=0, leaveout=1,
              interpolate = FALSE, confint = TRUE, alpha = 0.05,
              savefuns=FALSE, savepatterns=FALSE, verbose=FALSE)
}
\arguments{
  \item{X}{
    Either a point pattern (object of class \code{"ppp"}, \code{"lpp"}
    or other class), a fitted point process model (object of class \code{"ppm"},
    \code{"kppm"} or other class) or an envelope object (class
    \code{"envelope"}). 
  }
  \item{fun}{
    Function that computes the desired summary statistic
    for a point pattern.
  }
  \item{\dots}{
    Arguments passed to \code{\link{envelope}}.
  }
  \item{exponent}{
    Positive number. Exponent used in the test statistic. Use \code{exponent=2}
    for the Diggle-Cressie-Loosmore-Ford test, and \code{exponent=Inf}
    for the Maximum Absolute Deviation test.
    See Details.
  }
  \item{nsim}{
    Number of repetitions of the basic test.
  }
  \item{nsimsub}{
    Number of simulations in each basic test. There will be \code{nsim}
    repetitions of the basic test, each involving \code{nsimsub} simulated
    realisations, so there will be a total
    of \code{nsim * (nsimsub + 1)} simulations.
  }
  \item{alternative}{
    Character string specifying the alternative hypothesis.
    The default (\code{alternative="two.sided"}) is that the
    true value of the summary function is not equal to the theoretical
    value postulated under the null hypothesis.
    If \code{alternative="less"} the alternative hypothesis is that the
    true value of the summary function is lower than the theoretical value.
  }
  \item{rmin}{
    Optional. Left endpoint for the interval of \eqn{r} values
    on which the test statistic is calculated.
  }
  \item{leaveout}{
    Optional integer 0, 1 or 2 indicating how to calculate the
    deviation between the observed summary function and the
    nominal reference value, when the reference value must be estimated
    by simulation. See Details.
  }
  \item{interpolate}{
    Logical value indicating whether to interpolate the distribution of
    the test statistic by kernel smoothing, as described in
    Dao and Genton (2014, Section 5).
  }
  \item{confint}{
    Logical value indicating whether to compute a confidence interval
    for the \sQuote{true} \eqn{p}-value.
  }
  \item{alpha}{
    Significance level to be plotted (this has no effect on the calculation
    but is simply plotted as a reference value).
  }
  \item{savefuns}{
    Logical flag indicating whether to save the simulated
    function values (from the first stage).
  }
  \item{savepatterns}{
    Logical flag indicating whether to save the simulated
    point patterns (from the first stage). 
  }
  \item{verbose}{
    Logical flag indicating whether to print progress reports.
  }
}
\details{
  The Dao and Genton (2014) test for a spatial point pattern
  is described in \code{\link{dg.test}}.
  This test depends on the choice of an interval of
  distance values (the argument \code{rinterval}).
  A \emph{significance trace} (Bowman and Azzalini, 1997;
  Baddeley et al, 2014, 2015)
  of the test is a plot of the \eqn{p}-value
  obtained from the test against the length of
  the interval \code{rinterval}.
  
  The command \code{dg.sigtrace} effectively performs 
  \code{\link{dg.test}} on \code{X} using all possible intervals
  of the form \eqn{[0,R]}, and returns the resulting \eqn{p}-values
  as a function of \eqn{R}.

  The result is an object of class \code{"fv"} that can be plotted to
  obtain the significance trace. The plot shows the
  Dao-Genton adjusted
  \eqn{p}-value (solid black line), 
  the critical value \code{0.05} (dashed red line),
  and a pointwise 95\% confidence band (grey shading)
  for the \sQuote{true} (Neyman-Pearson) \eqn{p}-value.
  The confidence band is based on the Agresti-Coull (1998)
  confidence interval for a binomial proportion.

  If \code{X} is an envelope object and \code{fun=NULL} then
  the code will re-use the simulated functions stored in \code{X}.

  If the argument \code{rmin} is given, it specifies the left endpoint
  of the interval defining the test statistic: the tests are
  performed using intervals \eqn{[r_{\mbox{\scriptsize min}},R]}{[rmin,R]}
  where \eqn{R \ge r_{\mbox{\scriptsize min}}}{R >= rmin}.

  The argument \code{leaveout} specifies how to calculate the
  discrepancy between the summary function for the data and the
  nominal reference value, when the reference value must be estimated
  by simulation. The values \code{leaveout=0} and
  \code{leaveout=1} are both algebraically equivalent (Baddeley et al, 2014,
  Appendix) to computing the difference \code{observed - reference}
  where the \code{reference} is the mean of simulated values.
  The value \code{leaveout=2} gives the leave-two-out discrepancy
  proposed by Dao and Genton (2014).
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

  Dao, N.A. and Genton, M. (2014)
  A Monte Carlo adjusted goodness-of-fit test for
  parametric models describing spatial point patterns.
  \emph{Journal of Graphical and Computational Statistics}
  \bold{23}, 497--517.
}
\author{
  \adrian, Andrew Hardegen, Tom Lawrence,
  Robin Milne, Gopalan Nair and Suman Rakshit.
  Implemented by
  \spatstatAuthors.
}
\seealso{
  \code{\link{dg.test}} for the Dao-Genton test,
  \code{\link{dclf.sigtrace}} for significance traces of other tests.
}
\examples{
  ns <- if(interactive()) 19 else 5
  plot(dg.sigtrace(cells, nsim=ns))
}
\keyword{spatial}
\keyword{htest}
