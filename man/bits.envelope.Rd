\name{bits.envelope}
\alias{bits.envelope}
\title{
  Global Envelopes for Balanced Independent Two-Stage Test
}
\description{
  Computes the global envelopes
  corresponding to the balanced independent two-stage
  Monte Carlo test of goodness-of-fit.
}
\usage{
bits.envelope(X, \dots,
            nsim = 19, nrank = 1,
            alternative=c("two.sided", "less", "greater"),
            leaveout=1, interpolate = FALSE,
            savefuns=FALSE, savepatterns=FALSE,
            verbose = TRUE)
}
\arguments{
  \item{X}{
    Either a point pattern dataset (object of class \code{"ppp"},
    \code{"lpp"} or \code{"pp3"}) or a fitted point process model
    (object of class \code{"ppm"}, \code{"kppm"} or \code{"slrm"}).
  }
  \item{\dots}{
    Arguments passed to 
    \code{\link{mad.test}} or \code{\link{envelope}} to
    control the conduct of the test.
    Useful arguments include \code{fun} to determine the summary
    function, \code{rinterval} to determine the range of
    \eqn{r} values used in the test, and
    \code{verbose=FALSE} to turn off the messages.
  }
  \item{nsim}{
    Number of simulated patterns to be generated in each stage.
    Number of simulations in each basic test. There will be \code{nsim}
    repetitions of the basic test, each involving \code{nsim} simulated
    realisations, together with one independent set of \code{nsim}
    realisations, so there will be a total
    of \code{nsim * (nsim + 1)} simulations.
  }
  \item{nrank}{
    Integer. Rank of the envelope value amongst the \code{nsim} simulated
    values. A rank of 1 means that the minimum and maximum
    simulated values will be used.
  }
  \item{alternative}{
    Character string determining whether the envelope corresponds
    to a two-sided test (\code{alternative="two.sided"}, the default)
    or a one-sided test with a lower critical boundary
    (\code{alternative="less"}) or a one-sided test
    with an upper critical boundary (\code{alternative="greater"}).
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
    Logical value determining whether to print progress reports.
  }
}
\details{
  Computes global simulation envelopes corresponding to the
  balanced independent two-stage Monte Carlo test of goodness-of-fit
  described by Baddeley et al (2017). 
  The envelopes are described in Baddeley et al (2019).

  If \code{X} is a point pattern, the null hypothesis is CSR.

  If \code{X} is a fitted model, the null hypothesis is that model.

  This command is similar to \code{\link{dg.envelope}} which corresponds
  to the Dao-Genton test of goodness-of-fit.
  It was shown in Baddeley et al (2017) that
  the Dao-Genton test is biased when the significance level is very small
  (small \eqn{p}-values are not reliable) and
  we recommend \code{\link{bits.envelope}} in this case.
}
\value{
  An object of class \code{"fv"}.
}
\references{
  Dao, N.A. and Genton, M. (2014)
  A Monte Carlo adjusted goodness-of-fit test for
  parametric models describing spatial point patterns.
  \emph{Journal of Graphical and Computational Statistics}
  \bold{23}, 497--517.

  Baddeley, A., Hardegen, A., Lawrence, T.,
  Milne, R.K., Nair, G. and Rakshit, S. (2017)
  On two-stage Monte Carlo tests of composite
  hypotheses. \emph{Computational Statistics and Data Analysis}
  \bold{114}, {75--87}.

  Baddeley, A., Hardegen, A., Lawrence, L., 
  Milne, R.K., Nair, G.M. and Rakshit, S. (2019)
  Pushing the envelope: extensions of graphical
  Monte Carlo tests. In preparation.
}
\author{
  Adrian Baddeley, Andrew Hardegen, Tom Lawrence,
  Robin Milne, Gopalan Nair and Suman Rakshit.
  Implemented by \adrian.
}
\seealso{
  \code{\link{dg.envelope}},
  \code{\link{bits.test}},
  \code{\link{mad.test}},
  \code{\link{envelope}}
}
\examples{
  ns <- if(interactive()) 19 else 4
  E <- bits.envelope(swedishpines, Lest, nsim=ns)
  E
  plot(E)
  Eo <- bits.envelope(swedishpines, Lest, alternative="less", nsim=ns)
  Ei <- bits.envelope(swedishpines, Lest, interpolate=TRUE, nsim=ns)
}
\keyword{spatial}
\keyword{htest}
\keyword{hplot}
\keyword{iteration}
