\name{bits.test}
\alias{bits.test}
\title{
  Balanced Independent Two-Stage Monte Carlo Test
}
\description{
  Performs a Balanced Independent Two-Stage Monte Carlo test
  of goodness-of-fit for spatial pattern.
}
\usage{
bits.test(X, \dots,
        exponent = 2, nsim=19, 
        alternative=c("two.sided", "less", "greater"),
        leaveout=1, interpolate = FALSE,
        savefuns=FALSE, savepatterns=FALSE,
        verbose = TRUE)
}
\arguments{
  \item{X}{
    Either a point pattern dataset (object of class \code{"ppp"},
    \code{"lpp"} or \code{"pp3"}) or a fitted point process model
    (object of class \code{"ppm"}, \code{"kppm"}, \code{"lppm"}
    or \code{"slrm"}).
  }
  \item{\dots}{
    Arguments passed to \code{\link{dclf.test}} or
    \code{\link{mad.test}} or \code{\link{envelope}} to
    control the conduct of the test.
    Useful arguments include \code{fun} to determine the summary
    function, \code{rinterval} to determine the range of
    \eqn{r} values used in the test, 
    and \code{use.theory} described under Details. 
  }
  \item{exponent}{
    Exponent used in the test statistic. Use \code{exponent=2}
    for the Diggle-Cressie-Loosmore-Ford test, and \code{exponent=Inf}
    for the Maximum Absolute Deviation test.
  }
  \item{nsim}{
    Number of replicates in each stage of the test.
    A total of \code{nsim * (nsim + 1)} simulated point patterns will be
    generated, and the \eqn{p}-value will be a multiple of \code{1/(nsim+1)}.
  }
  \item{alternative}{
    Character string specifying the alternative hypothesis.
    The default (\code{alternative="two.sided"}) is that the
    true value of the summary function is not equal to the theoretical
    value postulated under the null hypothesis.
    If \code{alternative="less"} the alternative hypothesis is that the
    true value of the summary function is lower than the theoretical value.
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
  \item{savefuns}{
    Logical flag indicating whether to save the simulated
    function values (from the first stage).
  }
  \item{savepatterns}{
    Logical flag indicating whether to save the simulated
    point patterns (from the first stage). 
  }
  \item{verbose}{
    Logical value indicating whether to print progress reports.
  }
}
\details{
  Performs the Balanced Independent Two-Stage Monte Carlo test
  proposed by Baddeley et al (2017),
  an improvement of the Dao-Genton (2014) test.
  
  If \code{X} is a point pattern, the null hypothesis is CSR.

  If \code{X} is a fitted model, the null hypothesis is that model.

  The argument \code{use.theory} passed to \code{\link{envelope}}
  determines whether to compare the summary function for the data
  to its theoretical value for CSR (\code{use.theory=TRUE})
  or to the sample mean of simulations from CSR
  (\code{use.theory=FALSE}).

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
  A hypothesis test (object of class \code{"htest"}
  which can be printed to show the outcome of the test.
}
\references{
  Dao, N.A. and Genton, M. (2014)
  A Monte Carlo adjusted goodness-of-fit test for
  parametric models describing spatial point patterns.
  \emph{Journal of Graphical and Computational Statistics}
  \bold{23}, 497--517.

  Baddeley, A., Diggle, P.J., Hardegen, A., Lawrence, T., Milne,
  R.K. and Nair, G. (2014) On tests of spatial pattern based on
  simulation envelopes. \emph{Ecological Monographs} \bold{84} (3) 477--489.
  
  Baddeley, A., Hardegen, A., Lawrence, L., 
  Milne, R.K., Nair, G.M. and Rakshit, S. (2017)
  On two-stage Monte Carlo tests of composite hypotheses.
  \emph{Computational Statistics and Data Analysis}, in press.
}
\author{
  Adrian Baddeley, Andrew Hardegen, Tom Lawrence,
  Robin Milne, Gopalan Nair and Suman Rakshit.
  Implemented by \spatstatAuthors.
}
\seealso{
  \code{\link{dg.test}},
  \code{\link{dclf.test}},
  \code{\link{mad.test}}
}
\examples{
 ns <- if(interactive()) 19 else 4
 bits.test(cells, nsim=ns)
 bits.test(cells, alternative="less", nsim=ns)
 bits.test(cells, nsim=ns, interpolate=TRUE)
}
\keyword{spatial}
\keyword{htest}
