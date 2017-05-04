\name{envelope}
\alias{envelope}
\alias{envelope.ppp}
\alias{envelope.ppm}
\alias{envelope.kppm}
\title{Simulation Envelopes of Summary Function}
\description{
  Computes simulation envelopes of a summary function.
}
\usage{
  envelope(Y, fun, \dots)

  \method{envelope}{ppp}(Y, fun=Kest, nsim=99, nrank=1, \dots,
  funargs=list(), funYargs=funargs,
  simulate=NULL,  fix.n=FALSE, fix.marks=FALSE,
  verbose=TRUE, clipdata=TRUE,
  transform=NULL, global=FALSE, ginterval=NULL, use.theory=NULL, 
  alternative=c("two.sided", "less", "greater"),
  scale=NULL, clamp=FALSE, 
  savefuns=FALSE, savepatterns=FALSE,
  nsim2=nsim, VARIANCE=FALSE, nSD=2, Yname=NULL, maxnerr=nsim,
  do.pwrong=FALSE, envir.simul=NULL)

  \method{envelope}{ppm}(Y, fun=Kest, nsim=99, nrank=1, \dots, 
  funargs=list(), funYargs=funargs,
  simulate=NULL, fix.n=FALSE, fix.marks=FALSE,
  verbose=TRUE, clipdata=TRUE,
  start=NULL, control=update(default.rmhcontrol(Y), nrep=nrep), nrep=1e5,
  transform=NULL, global=FALSE, ginterval=NULL, use.theory=NULL, 
  alternative=c("two.sided", "less", "greater"), 
  scale=NULL, clamp=FALSE, 
  savefuns=FALSE, savepatterns=FALSE,
  nsim2=nsim, VARIANCE=FALSE, nSD=2, Yname=NULL, maxnerr=nsim,
  do.pwrong=FALSE, envir.simul=NULL)

  \method{envelope}{kppm}(Y, fun=Kest, nsim=99, nrank=1, \dots, 
  funargs=list(), funYargs=funargs,
  simulate=NULL,
  verbose=TRUE, clipdata=TRUE,
  transform=NULL, global=FALSE, ginterval=NULL, use.theory=NULL, 
  alternative=c("two.sided", "less", "greater"), 
  scale=NULL, clamp=FALSE, 
  savefuns=FALSE, savepatterns=FALSE,
  nsim2=nsim, VARIANCE=FALSE, nSD=2, Yname=NULL, maxnerr=nsim,
  do.pwrong=FALSE, envir.simul=NULL)
}
\arguments{
  \item{Y}{
    Object containing point pattern data.
    A point pattern (object of class
    \code{"ppp"}) or a fitted point process model
    (object of class \code{"ppm"} or \code{"kppm"}).
  }
  \item{fun}{
    Function that computes the desired summary statistic
    for a point pattern. 
  }
  \item{nsim}{
    Number of simulated point patterns to be generated
    when computing the envelopes.
  }
  \item{nrank}{
    Integer. Rank of the envelope value amongst the \code{nsim} simulated
    values. A rank of 1 means that the minimum and maximum
    simulated values will be used.
  }
  \item{\dots}{
    Extra arguments passed to \code{fun}.
  }
  \item{funargs}{
    A list, containing extra arguments to be passed to \code{fun}.
  }
  \item{funYargs}{
    Optional. A list, containing extra arguments to be passed to
    \code{fun} when applied to the original data \code{Y} only.
  }
  \item{simulate}{
    Optional. Specifies how to generate the simulated point patterns.
    If \code{simulate} is an expression in the R language, then this
    expression will be evaluated \code{nsim} times,
    to obtain \code{nsim} point patterns which are taken as the
    simulated patterns from which the envelopes are computed.
    If \code{simulate} is a list of point patterns, then the entries
    in this list will be treated as the simulated patterns from which
    the envelopes are computed.
    Alternatively \code{simulate} may be an object produced by the
    \code{envelope} command: see Details.
  }
  \item{fix.n}{
    Logical. If \code{TRUE}, simulated patterns will have the
    same number of points as the original data pattern.
    This option is currently not available for \code{envelope.kppm}.
  }
  \item{fix.marks}{
    Logical. If \code{TRUE}, simulated patterns will have the
    same number of points \emph{and} the same marks as the
    original data pattern. In a multitype point pattern this means that
    the simulated patterns will have the same number of points
    \emph{of each type} as the original data.
    This option is currently not available for \code{envelope.kppm}.
  }
  \item{verbose}{
    Logical flag indicating whether to print progress reports
    during the simulations.
  }
  \item{clipdata}{
    Logical flag indicating whether the data point pattern should be
    clipped to the same window as the simulated patterns,
    before the summary function for the data is computed.
    This should usually be \code{TRUE} to ensure that the
    data and simulations are properly comparable.
  }
  \item{start,control}{
    Optional. These specify the arguments \code{start} and \code{control}
    of \code{rmh}, giving complete control over the simulation
    algorithm. Applicable only when \code{Y} is a fitted model
    of class \code{"ppm"}.
  }
  \item{nrep}{
    Number of iterations in the Metropolis-Hastings simulation
    algorithm. Applicable only when \code{Y} is a fitted model
    of class \code{"ppm"}.
  }
  \item{transform}{
    Optional. A transformation to be applied to the
    function values, before the envelopes are computed.
    An expression object (see Details).
  }
  \item{global}{
    Logical flag indicating whether envelopes should be pointwise
    (\code{global=FALSE}) or simultaneous (\code{global=TRUE}).
  }
  \item{ginterval}{
    Optional.
    A vector of length 2 specifying
    the interval of \eqn{r} values for the simultaneous critical
    envelopes. Only relevant if \code{global=TRUE}.
  }
  \item{use.theory}{
    Logical value indicating whether to use the theoretical value,
    computed by \code{fun}, as the reference value for simultaneous
    envelopes. Applicable only when \code{global=TRUE}.
    Default is \code{use.theory=TRUE} if \code{Y} is a point pattern,
    or a point process model equivalent to Complete Spatial Randomness,
    and \code{use.theory=FALSE} otherwise.
  }
  \item{alternative}{
    Character string determining whether the envelope corresponds
    to a two-sided test (\code{side="two.sided"}, the default)
    or a one-sided test with a lower critical boundary
    (\code{side="less"}) or a one-sided test
    with an upper critical boundary (\code{side="greater"}).
  }
  \item{scale}{
    Optional. Scaling function for global envelopes.
    A function in the \R language which determines the
    relative scale of deviations, as a function of
    distance \eqn{r}, when computing the global envelopes.
    Applicable only when \code{global=TRUE}.
    Summary function values for distance \code{r}
    will be \emph{divided} by \code{scale(r)} before the
    maximum deviation is computed. The resulting global envelopes
    will have width proportional to \code{scale(r)}. 
  }
  \item{clamp}{
    Logical value indicating how to compute envelopes when
    \code{alternative="less"} or \code{alternative="greater"}.
    Deviations of the observed
    summary function from the theoretical summary function are initially
    evaluated as signed real numbers, with large positive values indicating
    consistency with the alternative hypothesis.
    If \code{clamp=FALSE} (the default), these values are not changed.
    If \code{clamp=TRUE}, any negative values are replaced by zero.
  }
  \item{savefuns}{
    Logical flag indicating whether to save all the simulated
    function values.
  }
  \item{savepatterns}{
    Logical flag indicating whether to save all the simulated
    point patterns.
  }
  \item{nsim2}{
    Number of extra simulated point patterns to be generated
    if it is necessary to use simulation to estimate the theoretical
    mean of the summary function. Only relevant when \code{global=TRUE}
    and the simulations are not based on CSR.
  }
  \item{VARIANCE}{
    Logical. If \code{TRUE}, critical envelopes will be calculated
    as sample mean plus or minus \code{nSD} times sample standard
    deviation.
  }
  \item{nSD}{
    Number of estimated standard deviations used to determine
    the critical envelopes, if \code{VARIANCE=TRUE}.
  }
  \item{Yname}{
    Character string that should be used as the name of the 
    data point pattern \code{Y} when printing or plotting the results.
  }
  \item{maxnerr}{
    Maximum number of rejected patterns.
    If \code{fun} yields an error when applied to a simulated point
    pattern (for example, because the pattern is empty and \code{fun}
    requires at least one point), the pattern will be rejected
    and a new random point pattern will be generated. If this happens
    more than \code{maxnerr} times, the algorithm will give up.
  }
  \item{do.pwrong}{
    Logical. If \code{TRUE}, the algorithm will also estimate
    the true significance level of the \dQuote{wrong} test (the test that
    declares the summary function for the data to be significant
    if it lies outside the \emph{pointwise} critical boundary at any
    point). This estimate is printed when the result is printed.
  }
  \item{envir.simul}{
    Environment in which to evaluate the expression \code{simulate},
    if not the current environment.
  }
}
\value{
  An object of class \code{"envelope"}
  and \code{"fv"}, see \code{\link{fv.object}},
  which can be printed and plotted directly.

  Essentially a data frame containing columns
  \item{r}{the vector of values of the argument \eqn{r} 
    at which the summary function \code{fun} has been  estimated
  }
  \item{obs}{
    values of the summary function for the data point pattern
  }
  \item{lo}{
    lower envelope of simulations
  }
  \item{hi}{
    upper envelope of simulations
  }
  and \emph{either}
  \item{theo}{
    theoretical value of the summary function under CSR
    (Complete Spatial Randomness, a uniform Poisson point process)
    if the simulations were generated according to CSR
  }
  \item{mmean}{
    estimated theoretical value of the summary function,
    computed by averaging simulated values, 
    if the simulations were not generated according to CSR.
  }
  Additionally, if \code{savepatterns=TRUE}, the return value has an attribute
  \code{"simpatterns"} which is a list containing the \code{nsim}
  simulated patterns. If \code{savefuns=TRUE}, the return value
  has an attribute \code{"simfuns"} which is an object of class
  \code{"fv"} containing the summary functions
  computed for each of the \code{nsim} simulated patterns.
}
\details{
  The \code{envelope} command performs simulations and
  computes envelopes of a summary statistic based on the simulations.
  The result is an object that can be plotted to display the envelopes.
  The envelopes can be used to assess the goodness-of-fit of
  a point process model to point pattern data.

  For the most basic use, if you have a point pattern \code{X} and
  you want to test Complete Spatial Randomness (CSR), type
  \code{plot(envelope(X, Kest,nsim=39))} to see the \eqn{K} function
  for \code{X} plotted together with the envelopes of the
  \eqn{K} function for 39 simulations of CSR. 
  
  The \code{envelope} function is generic, with methods for
  the classes \code{"ppp"}, \code{"ppm"} and \code{"kppm"}
  described here. There are also methods for the classes \code{"pp3"},
  \code{"lpp"} and \code{"lppm"} which are described separately
  under \code{\link{envelope.pp3}} and \code{\link{envelope.lpp}}.
  Envelopes can also be computed from other envelopes, using
  \code{\link{envelope.envelope}}.
  
  To create simulation envelopes, the command \code{envelope(Y, ...)} 
  first generates \code{nsim} random point patterns
  in one of the following ways. 
  \itemize{
    \item 
    If \code{Y} is a point pattern (an object of class \code{"ppp"})
    and \code{simulate=NULL},
    then we generate \code{nsim} simulations of
    Complete Spatial Randomness (i.e. \code{nsim} simulated point patterns
    each being a realisation of the uniform Poisson point process)
    with the same intensity as the pattern \code{Y}.
    (If \code{Y} is a multitype point pattern, then the simulated patterns
    are also given independent random marks; the probability
    distribution of the random marks is determined by the
    relative frequencies of marks in \code{Y}.)
    \item
    If \code{Y} is a fitted point process model (an object of class
    \code{"ppm"} or \code{"kppm"}) and \code{simulate=NULL},
    then this routine generates \code{nsim} simulated
    realisations of that model.
    \item
    If \code{simulate} is supplied, then it determines how the
    simulated point patterns are generated. It may be either
    \itemize{
      \item
      an expression in the R language, typically containing a call
      to a random generator. This expression will be evaluated
      \code{nsim} times to yield \code{nsim} point patterns. For example
      if \code{simulate=expression(runifpoint(100))} then each simulated
      pattern consists of exactly 100 independent uniform random points.
      \item
      a list of point patterns.
      The entries in this list will be taken as the simulated patterns.
      \item
      an object of class \code{"envelope"}. This should have been
      produced by calling \code{envelope} with the
      argument \code{savepatterns=TRUE}.
      The simulated point patterns that were saved in this object
      will be extracted and used as the simulated patterns for the
      new envelope computation. This makes it possible to plot envelopes
      for two different summary functions based on exactly the same set of
      simulated point patterns.
      }
  }
  
  The summary statistic \code{fun} is applied to each of these simulated
  patterns. Typically \code{fun} is one of the functions
  \code{Kest}, \code{Gest}, \code{Fest}, \code{Jest}, \code{pcf},
  \code{Kcross}, \code{Kdot}, \code{Gcross}, \code{Gdot},
  \code{Jcross}, \code{Jdot}, \code{Kmulti}, \code{Gmulti},
  \code{Jmulti} or \code{Kinhom}. It may also be a character string
  containing the name of one of these functions.

  The statistic \code{fun} can also be a user-supplied function;
  if so, then it must have arguments \code{X} and \code{r}
  like those in the functions listed above, and it must return an object
  of class \code{"fv"}.

  Upper and lower critical envelopes are computed in one of the following ways:
  \describe{
    \item{pointwise:}{by default, envelopes are calculated pointwise
      (i.e. for each value of the distance argument \eqn{r}), by sorting the
      \code{nsim} simulated values, and taking the \code{m}-th lowest
      and \code{m}-th highest values, where \code{m = nrank}.
      For example if \code{nrank=1}, the upper and lower envelopes
      are the pointwise maximum and minimum of the simulated values.

      The pointwise envelopes are \bold{not} \dQuote{confidence bands}
      for the true value of the function! Rather,
      they specify the critical points for a Monte Carlo test
      (Ripley, 1981). The test is constructed by choosing a
      \emph{fixed} value of \eqn{r}, and rejecting the null hypothesis if the
      observed function value
      lies outside the envelope \emph{at this value of} \eqn{r}.
      This test has exact significance level
      \code{alpha = 2 * nrank/(1 + nsim)}.
    }
    \item{simultaneous:}{if \code{global=TRUE}, then the envelopes are
      determined as follows. First we calculate the theoretical mean value of
      the summary statistic (if we are testing CSR, the theoretical
      value is supplied by \code{fun}; otherwise we perform a separate
      set of \code{nsim2} simulations, compute the
      average of all these simulated values, and take this average
      as an estimate of the theoretical mean value). Then, for each simulation,
      we compare the simulated curve to the theoretical curve, and compute the
      maximum absolute difference between them (over the interval
      of \eqn{r} values specified by \code{ginterval}). This gives a
      deviation value \eqn{d_i}{d[i]} for each of the \code{nsim}
      simulations. Finally we take the \code{m}-th largest of the
      deviation values, where \code{m=nrank}, and call this
      \code{dcrit}. Then the simultaneous envelopes are of the form
      \code{lo = expected - dcrit} and \code{hi = expected + dcrit} where
      \code{expected} is either the theoretical mean value \code{theo}
      (if we are testing CSR) or the estimated theoretical value
      \code{mmean} (if we are testing another model). The simultaneous critical
      envelopes have constant width \code{2 * dcrit}.

      The simultaneous critical envelopes allow us to perform a different
      Monte Carlo test (Ripley, 1981). The test rejects the null
      hypothesis if the graph of the observed function
      lies outside the envelope \bold{at any value of} \eqn{r}.
      This test has exact significance level
      \code{alpha = nrank/(1 + nsim)}.

      This test can also be performed using \code{\link{mad.test}}.
    }
    \item{based on sample moments:}{if \code{VARIANCE=TRUE},
      the algorithm calculates the
      (pointwise) sample mean and sample variance of
      the simulated functions. Then the envelopes are computed
      as mean plus or minus \code{nSD} standard deviations.
      These envelopes do not have an exact significance interpretation.
      They are a naive approximation to
      the critical points of the Neyman-Pearson test
      assuming the summary statistic is approximately Normally
      distributed.
    }
  }
  
  The return value is an object of class \code{"fv"} containing
  the summary function for the data point pattern,
  the upper and lower simulation envelopes, and 
  the theoretical expected value (exact or estimated) of the summary function 
  for the model being tested. It can be plotted
  using \code{\link{plot.envelope}}.

  If \code{VARIANCE=TRUE} then the return value also includes the
  sample mean, sample variance and other quantities.

  Arguments can be passed to the function \code{fun} through
  \code{...}. This means that you simply specify these arguments in the call to
  \code{envelope}, and they will be passed to \code{fun}.
  In particular, the argument \code{correction}
  determines the edge correction to be used to calculate the summary
  statistic. See the section on Edge Corrections, and the Examples.

  Arguments can also be passed to the function \code{fun}
  through the list \code{funargs}. This mechanism is typically used if
  an argument of \code{fun} has the same name as an argument of
  \code{envelope}. The list \code{funargs} should contain
  entries of the form \code{name=value}, where each \code{name} is the name
  of an argument of \code{fun}.

  There is also an option, rarely used, in which different function
  arguments are used when computing the summary function
  for the data \code{Y} and for the simulated patterns.
  If \code{funYargs} is given, it will be used
  when the summary function for the data \code{Y} is computed,
  while \code{funargs} will be used when computing the summary function
  for the simulated patterns.
  This option is only needed in rare cases: usually the basic principle
  requires that the data and simulated patterns must be treated
  equally, so that \code{funargs} and \code{funYargs} should be identical.

  If \code{Y} is a fitted cluster point process model (object of
  class \code{"kppm"}), and \code{simulate=NULL},
  then the model is simulated directly
  using \code{\link{simulate.kppm}}.
  
  If \code{Y} is a fitted Gibbs point process model (object of
  class \code{"ppm"}), and \code{simulate=NULL},
  then the model is simulated
  by running the Metropolis-Hastings algorithm \code{\link{rmh}}.
  Complete control over this algorithm is provided by the 
  arguments \code{start} and \code{control} which are passed
  to \code{\link{rmh}}.

  For simultaneous critical envelopes (\code{global=TRUE})
  the following options are also useful:
  \describe{
    \item{\code{ginterval}}{determines the interval of \eqn{r} values
      over which the deviation between curves is calculated.
      It should be a numeric vector of length 2.
      There is a sensible default (namely, the recommended plotting
      interval for \code{fun(X)}, or the range of \code{r} values if
      \code{r} is explicitly specified).
    }
    \item{\code{transform}}{specifies a transformation of the
      summary function \code{fun} that will be carried out before the
      deviations are computed.
      Such transforms are useful if \code{global=TRUE} or
      \code{VARIANCE=TRUE}.
      The \code{transform} must be an expression object
      using the symbol \code{.} to represent the function value
      (and possibly other symbols recognised by \code{\link{with.fv}}).
      For example, 
      the conventional way to normalise the \eqn{K} function
      (Ripley, 1981) is to transform it to the \eqn{L} function
      \eqn{L(r) = \sqrt{K(r)/\pi}}{L(r) = sqrt(K(r)/\pi)}
      and this is implemented by setting
      \code{transform=expression(sqrt(./pi))}.
    }
  }

  It is also possible to extract the summary functions for each of the
  individual simulated point patterns, by setting \code{savefuns=TRUE}.
  Then the return value also 
  has an attribute \code{"simfuns"} containing all the 
  summary functions for the individual simulated patterns.
  It is an \code{"fv"} object containing
  functions named \code{sim1, sim2, ...} representing the \code{nsim}
  summary functions.

  It is also possible to save the simulated point patterns themselves,
  by setting \code{savepatterns=TRUE}. Then the return value also has
  an attribute \code{"simpatterns"} which is a list of length
  \code{nsim} containing all the simulated point patterns.

  See \code{\link{plot.envelope}} and \code{\link{plot.fv}}
  for information about how to plot the envelopes.

  Different envelopes can be recomputed from the same data
  using \code{\link{envelope.envelope}}.
  Envelopes can be combined using \code{\link{pool.envelope}}.
}
\section{Errors and warnings}{
  An error may be generated if one of the simulations produces a
  point pattern that is empty, or is otherwise unacceptable to the
  function \code{fun}.
  
  The upper envelope may be \code{NA} (plotted as plus or minus
  infinity) if some of the function values
  computed for the simulated point patterns are \code{NA}.
  Whether this occurs will depend on the function \code{fun},
  but it usually happens when the simulated point pattern does not contain
  enough points to compute a meaningful value.
}
\section{Confidence intervals}{
  Simulation envelopes do \bold{not} compute confidence intervals;
  they generate significance bands. 
  If you really need a confidence interval for the true summary function
  of the point process, use \code{\link{lohboot}}.
  See also \code{\link{varblock}}.
}
\section{Edge corrections}{
  It is common to apply a correction for edge effects when
  calculating a summary function such as the \eqn{K} function.
  Typically the user has a choice between several possible edge
  corrections.
  In a call to \code{envelope}, the user can specify the edge correction
  to be applied in \code{fun}, using the argument \code{correction}.
  See the Examples below.

  \describe{
    \item{Summary functions in \pkg{spatstat}}{
      Summary functions that are available in \pkg{spatstat}, such as
      \code{\link{Kest}}, \code{\link{Gest}} and \code{\link{pcf}},
      have a standard argument called \code{correction} which specifies
      the name of one or more edge corrections.

      The list of available edge
      corrections is different for each summary function,
      and may also depend on the kind of window in which the point pattern is
      recorded.
      In the
      case of \code{Kest} (the default and most frequently used value of
      \code{fun}) the best edge correction is Ripley's isotropic
      correction if the window is rectangular or polygonal, 
      and the translation correction if the window is a binary mask.
      See the help files for the individual
      functions for more information.
      
      All the summary functions in \pkg{spatstat}
      recognise the option \code{correction="best"}
      which gives the \dQuote{best} (most accurate) available edge correction
      for that function. 

      In a call to \code{envelope}, if \code{fun} is one of the
      summary functions provided in \pkg{spatstat}, then the default
      is \code{correction="best"}. This means that
      \emph{by default, the envelope will be computed
	using the \dQuote{best} available edge correction}.

      The user can override this default by specifying the argument
      \code{correction}. For example the computation can be accelerated
      by choosing another edge correction which is less accurate
      than the \dQuote{best} one, but faster to compute.
    }
    \item{User-written summary functions}{
      If \code{fun} is a function written by the user,
      then \code{envelope} has to guess what to do.

      If \code{fun} has an argument
      called \code{correction}, or has \code{\dots} arguments,
      then \code{envelope} assumes that the function 
      can handle a correction argument. To compute the envelope,
      \code{fun} will be called with a \code{correction} argument.
      The default is \code{correction="best"}, unless
      overridden in the call to \code{envelope}.

      Otherwise, if \code{fun} does not have an argument
      called \code{correction} and does not have \code{\dots} arguments,
      then \code{envelope} assumes that the function 
      \emph{cannot} handle a correction argument. To compute the
      envelope, \code{fun} is called without a correction argument.
    }
  }
}
\references{
  Baddeley, A., Diggle, P.J., Hardegen, A., Lawrence, T., Milne,
  R.K. and Nair, G. (2014) On tests of spatial pattern based on
  simulation envelopes. \emph{Ecological Monographs} \bold{84} (3) 477--489.
  
  Cressie, N.A.C. \emph{Statistics for spatial data}.
    John Wiley and Sons, 1991.

  Diggle, P.J. \emph{Statistical analysis of spatial point patterns}.
  Arnold, 2003.

  Ripley, B.D. (1981) 
  \emph{Spatial statistics}.
  John Wiley and Sons.

  Ripley, B.D. \emph{Statistical inference for spatial processes}.
  Cambridge University Press, 1988.

  Stoyan, D. and Stoyan, H. (1994)
  Fractals, random shapes and point fields:
  methods of geometrical statistics.
  John Wiley and Sons.
} 
\seealso{
  \code{\link{dclf.test}},
  \code{\link{mad.test}}
  for envelope-based tests.
  
  \code{\link{fv.object}},
  \code{\link{plot.envelope}},
  \code{\link{plot.fv}},
  \code{\link{envelope.envelope}},
  \code{\link{pool.envelope}}
  for handling envelopes.
  There are also methods for \code{print} and \code{summary}.
  
  \code{\link{Kest}},
  \code{\link{Gest}},
  \code{\link{Fest}},
  \code{\link{Jest}},
  \code{\link{pcf}},
  \code{\link{ppp}},
  \code{\link{ppm}},
  \code{\link{default.expand}}
}
\examples{
 X <- simdat

 # Envelope of K function under CSR
 \dontrun{
 plot(envelope(X))
 }
 \testonly{
  plot(envelope(X, nsim=3))
 }

 # Translation edge correction (this is also FASTER):
 \dontrun{
 plot(envelope(X, correction="translate"))
 }
 \testonly{
  E <- envelope(X, nsim=3, correction="translate")
 }

# Global envelopes
 \dontrun{
 plot(envelope(X, Lest, global=TRUE))
 plot(envelope(X, Kest,  global=TRUE, scale=function(r) { r }))
 }
 \testonly{
  E <- envelope(X, Lest, nsim=3, global=TRUE)
  E <- envelope(X, Kest, nsim=3, global=TRUE, scale=function(r) { r })
  E
  summary(E)
 }

 # Envelope of K function for simulations from Gibbs model 
 \dontrun{
 fit <- ppm(cells ~1, Strauss(0.05))
 plot(envelope(fit))
 plot(envelope(fit), global=TRUE)
 }
 \testonly{
  fit <- ppm(cells ~1, Strauss(0.05), nd=20)
  E <- envelope(fit, nsim=3, correction="border", nrep=100)
  E <- envelope(fit, nsim=3, correction="border", global=TRUE, nrep=100)
 }

 # Envelope of K function for simulations from cluster model 
 fit <- kppm(redwood ~1, "Thomas")
 \dontrun{
 plot(envelope(fit, Gest))
 plot(envelope(fit, Gest, global=TRUE))
 }
 \testonly{
  E <- envelope(fit, Gest, correction="rs", nsim=3, global=TRUE, nrep=100)
 }

 # Envelope of G function under CSR
 \dontrun{
 plot(envelope(X, Gest))
 }
 \testonly{
  E <- envelope(X, Gest, correction="rs", nsim=3)
 }

 # Envelope of L function under CSR
 #  L(r) = sqrt(K(r)/pi)
 \dontrun{
  E <- envelope(X, Kest)
  plot(E, sqrt(./pi) ~ r)
 }
 \testonly{
  E <- envelope(X, Kest, correction="border", nsim=3)
  plot(E, sqrt(./pi) ~ r)
 }

 # Simultaneous critical envelope for L function
 # (alternatively, use Lest)
 \dontrun{
  plot(envelope(X, Kest, transform=expression(sqrt(./pi)), global=TRUE))
 }
 \testonly{
  E <- envelope(X, Kest, nsim=3, correction="border",
               transform=expression(sqrt(./pi)), global=TRUE)
 }

 ## One-sided envelope
 \dontrun{
  plot(envelope(X, Lest, alternative="less"))
 }
 \testonly{
  E <- envelope(X, Lest, nsim=3, alternative="less")
 }
 
 # How to pass arguments needed to compute the summary functions:
 # We want envelopes for Jcross(X, "A", "B") 
 # where "A" and "B" are types of points in the dataset 'demopat'

 data(demopat)
 \dontrun{
 plot(envelope(demopat, Jcross, i="A", j="B"))
 }
 \testonly{
 plot(envelope(demopat, Jcross, correction="rs", i="A", j="B", nsim=3))
 }
 
 # Use of `simulate'
 \dontrun{
 plot(envelope(cells, Gest, simulate=expression(runifpoint(42))))
 plot(envelope(cells, Gest, simulate=expression(rMaternI(100,0.02))))
 }
 \testonly{
  plot(envelope(cells, Gest, correction="rs", simulate=expression(runifpoint(42)), nsim=3))
    plot(envelope(cells, Gest, correction="rs", simulate=expression(rMaternI(100, 0.02)),
nsim=3, global=TRUE))
 }

 # Envelope under random toroidal shifts
 data(amacrine)
 \dontrun{
 plot(envelope(amacrine, Kcross, i="on", j="off",
               simulate=expression(rshift(amacrine, radius=0.25)))) 
 }

 # Envelope under random shifts with erosion
 \dontrun{
 plot(envelope(amacrine, Kcross, i="on", j="off",
              simulate=expression(rshift(amacrine, radius=0.1, edge="erode"))))
 }
  
 # Envelope of INHOMOGENEOUS K-function with fitted trend

 # The following is valid.
 # Setting lambda=fit means that the fitted model is re-fitted to
 # each simulated pattern to obtain the intensity estimates for Kinhom.
 # (lambda=NULL would also be valid)

 fit <- kppm(redwood ~1, clusters="MatClust")
 \dontrun{
    plot(envelope(fit, Kinhom, lambda=fit, nsim=19))
 }
 \testonly{
    envelope(fit, Kinhom, lambda=fit, nsim=3)
 }

 # Note that the principle of symmetry, essential to the validity of
 # simulation envelopes, requires that both the observed and
 # simulated patterns be subjected to the same method of intensity
 # estimation. In the following example it would be incorrect to set the
 # argument 'lambda=red.dens' in the envelope command, because this
 # would mean that the inhomogeneous K functions of the simulated
 # patterns would be computed using the intensity function estimated
 # from the original redwood data, violating the symmetry.  There is
 # still a concern about the fact that the simulations are generated
 # from a model that was fitted to the data; this is only a problem in
 # small datasets.

\dontrun{
 red.dens <- density(redwood, sigma=bw.diggle)
 plot(envelope(redwood, Kinhom, sigma=bw.diggle,
         simulate=expression(rpoispp(red.dens))))
 }

 # Precomputed list of point patterns
\dontrun{
 nX <- npoints(X)
 PatList <- list()
 for(i in 1:19) PatList[[i]] <- runifpoint(nX)
 E <- envelope(X, Kest, nsim=19, simulate=PatList)
}
\testonly{
 PatList <- list()
 for(i in 1:3) PatList[[i]] <- runifpoint(10)
 E <- envelope(X, Kest, nsim=3, simulate=PatList)
}

# re-using the same point patterns
\dontrun{
 EK <- envelope(X, Kest, savepatterns=TRUE)
 EG <- envelope(X, Gest, simulate=EK)
}
\testonly{
 EK <- envelope(X, Kest, nsim=3, savepatterns=TRUE)
 EG <- envelope(X, Gest, nsim=3, simulate=EK)
}
}
\author{
  \spatstatAuthors.
}
\keyword{spatial}
\keyword{htest}
\keyword{hplot}
\keyword{iteration}
 
 
