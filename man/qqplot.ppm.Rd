\name{qqplot.ppm}
\alias{qqplot.ppm}
\title{
  Q-Q Plot of Residuals from Fitted Point Process Model
}
\description{
  Given a point process model fitted to a point pattern,
  produce a Q-Q plot based on residuals from the model.
}
\usage{
  qqplot.ppm(fit, nsim=100, expr=NULL, \dots, type="raw",
             style="mean", fast=TRUE, verbose=TRUE, plot.it=TRUE,
             dimyx=NULL, nrep=if(fast) 5e4 else 1e5,
             control=update(default.rmhcontrol(fit), nrep=nrep),
             saveall=FALSE,
             monochrome=FALSE,
             limcol=if(monochrome) "black" else "red",
             maxerr=max(100, ceiling(nsim/10)),
             check=TRUE, repair=TRUE, envir.expr)
}
\arguments{
  \item{fit}{
    The fitted point process model, which is to be assessed
    using the Q-Q plot. An object of class \code{"ppm"}.
    Smoothed residuals obtained from this fitted model will provide the
    ``data'' quantiles for the Q-Q plot.
  }
  \item{nsim}{
    The number of simulations from the ``reference'' point process model.
  }
  \item{expr}{
    Determines the simulation mechanism
    which provides the ``theoretical'' quantiles for the
    Q-Q plot. See Details.
  }
  \item{\dots}{
    Arguments passed to \code{\link{diagnose.ppm}} influencing the
    computation of residuals.
  }
  \item{type}{
    String indicating the type of residuals or weights to be used.
    Current options are \code{"eem"}
    for the Stoyan-Grabarnik exponential energy weights,
    \code{"raw"} for the raw residuals,
    \code{"inverse"} for the inverse-lambda residuals,
    and \code{"pearson"} for the Pearson residuals.
    A partial match is adequate.
  }
  \item{style}{
    Character string controlling the type of Q-Q plot.
    Options are \code{"classical"} and \code{"mean"}.
    See Details.
  }
  \item{fast}{
    Logical flag controlling the speed and accuracy of computation.
    Use \code{fast=TRUE} for interactive use and \code{fast=FALSE}
    for publication standard plots. See Details.
  }
  \item{verbose}{
    Logical flag controlling whether the algorithm prints progress
    reports during long computations.
  }
  \item{plot.it}{
    Logical flag controlling whether the function produces a plot
    or simply returns a value (silently).
  }
  \item{dimyx}{
    Dimensions of the pixel grid on which the smoothed residual field
    will be calculated. A vector of two integers.
  }
  \item{nrep}{
    If \code{control} is absent, then \code{nrep} gives the
    number of iterations of the Metropolis-Hastings algorithm
    that should be used to generate one simulation of the fitted point process.
  }
  \item{control}{
    List of parameters controlling the Metropolis-Hastings algorithm
    \code{\link{rmh}} which generates each simulated realisation from
    the model (unless the model is Poisson).
    This list becomes the argument \code{control}
    of \code{\link{rmh.default}}. It overrides \code{nrep}.
  }
  \item{saveall}{
    Logical flag indicating whether to save all the intermediate
    calculations. 
  }
  \item{monochrome}{
    Logical flag indicating whether the plot should be
    in black and white (\code{monochrome=TRUE}), or in colour
    (\code{monochrome=FALSE}).
  }
  \item{limcol}{
    String. The colour to be used when plotting the 95\% limit
    curves.
  }
  \item{maxerr}{
    Maximum number of failures tolerated while generating
    simulated realisations. See Details.
  }
  \item{check}{
    Logical value indicating whether to check the internal format
    of \code{fit}. If there is any possibility that this object
    has been restored from a dump file, or has otherwise lost track of
    the environment where it was originally computed, set
    \code{check=TRUE}. 
  }
  \item{repair}{
    Logical value indicating whether to repair the internal format
    of \code{fit}, if it is found to be damaged. 
  }
  \item{envir.expr}{
    Optional. An environment in which the expression
    \code{expr} should be evaluated.
  }
}
\value{
  An object of class \code{"qqppm"} containing the information
  needed to reproduce the Q-Q plot.
  Entries \code{x} and \code{y} are numeric vectors containing 
  quantiles of the simulations and of the data, respectively.
}
\details{
  This function generates a Q-Q plot of the residuals from a
  fitted point process model. It is an addendum to the suite of
  diagnostic plots produced by the function \code{\link{diagnose.ppm}},
  kept separate because it is computationally intensive. The
  quantiles of the theoretical distribution are estimated by simulation.

  In classical statistics, a Q-Q plot of residuals is a useful
  diagnostic for checking the distributional assumptions. Analogously,
  in spatial statistics, a Q-Q plot of the (smoothed) residuals from a
  fitted point process model is a useful way
  to check the interpoint interaction part of the model
  (Baddeley et al, 2005). The systematic part of the model
  (spatial trend, covariate effects, etc) is assessed using
  other plots made by \code{\link{diagnose.ppm}}.

  The argument \code{fit} represents the fitted point process
  model. It must be an object of class \code{"ppm"} (typically produced
  by the maximum pseudolikelihood fitting algorithm \code{\link{ppm}}).
  Residuals will be computed for this fitted model using
  \code{\link{residuals.ppm}},
  and the residuals will be kernel-smoothed to produce a ``residual
  field''. The values of this residual field will provide the
  ``data'' quantiles for the Q-Q plot.

  The argument \code{expr} is not usually specified.
  It provides a way to modify the ``theoretical'' or ``reference''
  quantiles for the Q-Q plot.

  In normal usage we set \code{expr=NULL}. The default
  is to generate \code{nsim} simulated realisations
  of the fitted model \code{fit}, re-fit this model to
  each of the simulated patterns,
  evaluate the residuals from
  these fitted models, and use the kernel-smoothed residual field
  from these fitted models as a sample from the reference distribution
  for the Q-Q plot.

  In advanced use, \code{expr} may be an \code{expression}.
  It will be re-evaluated \code{nsim} times, and should include
  random computations so that the results are not identical
  each time. The result of evaluating \code{expr}
  should be either a point pattern (object of class
  \code{"ppp"}) or a fitted point process model (object of class
  \code{"ppm"}). If the value is a point pattern, then the
  original fitted model \code{fit} will be fitted to this new point
  pattern using \code{\link{update.ppm}}, to yield another fitted
  model. Smoothed residuals obtained from these
  \code{nsim} fitted models will yield the ``theoretical'' quantiles for the
  Q-Q plot.

  Alternatively \code{expr} can be a list of point patterns,
  or an \code{envelope} object that contains a list of point patterns
  (typically generated by calling \code{\link{envelope}} with
  \code{savepatterns=TRUE}). These point patterns will be used
  as the simulated patterns. 

  Simulation is performed (if \code{expr=NULL})
  using the Metropolis-Hastings algorithm \code{\link{rmh}}.
  Each simulated realisation is the result of
  running the Metropolis-Hastings algorithm
  from an independent random starting state each time.
  The iterative and termination behaviour of the Metropolis-Hastings
  algorithm are governed by the argument \code{control}.
  See \code{\link{rmhcontrol}} for information about this argument.
  As a shortcut, the argument \code{nrep} determines
  the number of Metropolis-Hastings iterations used to generate
  each simulated realisation, if \code{control} is absent.

  By default, simulations are generated in an expanded
  window. Use the argument \code{control} to change this,
  as explained in the section on \emph{Warning messages}.
  
  The argument \code{type} selects the type of residual or weight
  that will be computed. For options, see \code{\link{diagnose.ppm}}.

  The argument \code{style} determines the type of Q-Q plot. It is
  highly recommended to use the default, \code{style="mean"}.
  \describe{
    \item{\code{style="classical"}}{
      The quantiles of the residual field for the data (on the \eqn{y}
      axis) are plotted against the
      quantiles of the \bold{pooled} simulations (on the \eqn{x} axis).
      This plot is biased, and therefore difficult to interpret,
      because of strong autocorrelations in the residual field
      and the large differences in sample size.
    }
    \item{\code{style="mean"}}{
      The order statistics of the residual field for the data are plotted
      against the sample means, over the \code{nsim} simulations,
      of the corresponding order statistics of the residual field
      for the simulated datasets.
      Dotted lines show the 2.5 and 97.5 percentiles, over the
      \code{nsim} simulations, of each order statistic.
    }
  }

  The argument \code{fast} is a simple way to control
  the accuracy and speed of computation.
  If \code{fast=FALSE}, the residual field is computed on
  a fine grid of pixels (by default 100 by 100 pixels, see below)
  and the Q-Q plot is based on the complete set of order statistics
  (usually 10,000 quantiles).
  If \code{fast=TRUE}, the residual field is computed on a coarse
  grid (at most 40 by 40 pixels) and the Q-Q plot is based on the
  \emph{percentiles} only. This is about 7 times faster.
  It is recommended to use \code{fast=TRUE} for interactive data
  analysis and \code{fast=FALSE} for definitive plots for
  publication.

  The argument \code{dimyx} gives full control over the resolution of the
  pixel grid used to calculate the smoothed residuals.
  Its interpretation is the same as the argument \code{dimyx}
  to the function \code{\link{as.mask}}.
  Note that \code{dimyx[1]} is the number of
  pixels in the \eqn{y} direction, and \code{dimyx[2]} is the number
  in the \eqn{x} direction. 
  If \code{dimyx} is not present, then the default pixel grid dimensions
  are controlled by \code{spatstat.options("npixel")}.

  Since the computation is so time-consuming, \code{qqplot.ppm} returns
  a list containing all the data necessary to re-display the Q-Q plot.
  It is advisable to assign the result of \code{qqplot.ppm} to something
  (or use \code{.Last.value} if you forgot to.)
  The return value is an object of class \code{"qqppm"}. There are methods for
  \code{\link{plot.qqppm}} and \code{\link{print.qqppm}}. See the
  Examples.

  The argument \code{saveall} is usually set to \code{FALSE}.
  If \code{saveall=TRUE}, then the intermediate results of calculation for each
  simulated realisation are saved and returned. The return value
  includes a 3-dimensional array \code{sim} containing the
  smoothed residual field images for each of the \code{nsim} 
  realisations. When \code{saveall=TRUE}, the return value is an object of very
  large size, and should not be saved on disk.

  Errors may occur during the simulation process, because
  random data are generated. For example:
  \itemize{
    \item one of the simulated patterns may be empty.
    \item one of the simulated patterns may 
    cause an error in the code that fits the point process model.
    \item the user-supplied argument \code{expr} may have a bug.
  }
  Empty point patterns do not cause a problem for the code,
  but they are reported.
  Other problems that would lead to a crash are trapped; 
  the offending simulated data are discarded, and the simulation is
  retried. The argument \code{maxerr} determines the maximum number of
  times that such errors will be tolerated (mainly as a
  safeguard against an infinite loop).
}

\section{Side Effects}{
	Produces a Q-Q plot if \code{plot.it} is TRUE.
}
\section{Warning messages}{
  A warning message will be issued if any of the simulations
  trapped an error (a potential crash).
  
  A warning message will be issued if all, or many, of the
  simulated point patterns are empty.
  This usually indicates a problem with the simulation procedure.
  
  The default behaviour of \code{qqplot.ppm} is to simulate patterns 
  on an expanded window (specified through the argument
  \code{control}) in order to avoid edge effects.
  The model's trend is extrapolated over this expanded
  window. If the trend is strongly inhomogeneous, the 
  extrapolated trend may have very large (or even infinite)
  values. This can cause the simulation algorithm to 
  produce empty patterns.

  The only way to suppress this problem entirely is to
  prohibit the expansion of the window, by setting
  the \code{control} argument to something like
  \code{control=list(nrep=1e6, expand=1)}. Here \code{expand=1}
  means there will be no expansion. See \code{\link{rmhcontrol}}
  for more information about the argument \code{control}.
}

\references{
  Baddeley, A., Turner, R., \ifelse{latex}{\out{M\o ller}}{Moller}, J. and Hazelton, M. (2005)
  Residual analysis for spatial point processes.
  \emph{Journal of the Royal Statistical Society, Series B}
  \bold{67}, 617--666.

  Stoyan, D. and Grabarnik, P. (1991)
  Second-order characteristics for stochastic structures connected with
  Gibbs point processes.
  \emph{Mathematische Nachrichten}, 151:95--100.
}
\seealso{
 \code{\link{diagnose.ppm}},
 \code{\link{lurking}},
 \code{\link{residuals.ppm}},
 \code{\link{eem}},
 \code{\link{ppm.object}},
 \code{\link{ppm}},
 \code{\link{rmh}},
 \code{\link{rmhcontrol}}
}
\examples{
    data(cells)

    fit <- ppm(cells, ~1, Poisson())
    diagnose.ppm(fit)  # no suggestion of departure from stationarity
    \dontrun{qqplot.ppm(fit, 80)  # strong evidence of non-Poisson interaction}
    \testonly{qqplot.ppm(fit, 4)}

    \dontrun{
     diagnose.ppm(fit, type="pearson")  
     qqplot.ppm(fit, type="pearson")
    }
    \testonly{qqplot.ppm(fit, 4, type="pearson")}

    ###########################################
    ## oops, I need the plot coordinates
    mypreciousdata <- .Last.value
    \dontrun{mypreciousdata <- qqplot.ppm(fit, type="pearson")}
    \testonly{mypreciousdata <- qqplot.ppm(fit, 4, type="pearson")}
    plot(mypreciousdata)

    ######################################################
    # Q-Q plots based on fixed n
    # The above QQ plots used simulations from the (fitted) Poisson process.
    # But I want to simulate conditional on n, instead of Poisson
    # Do this by setting rmhcontrol(p=1)
    fixit <- list(p=1)
    \dontrun{qqplot.ppm(fit, 100, control=fixit)}
    \testonly{qqplot.ppm(fit, 4, control=fixit)}

    ######################################################
    # Inhomogeneous Poisson data
    X <- rpoispp(function(x,y){1000 * exp(-3*x)}, 1000)
    plot(X)
    # Inhomogeneous Poisson model
    fit <- ppm(X, ~x, Poisson())
    \dontrun{qqplot.ppm(fit, 100)}
    \testonly{qqplot.ppm(fit, 4)}
    # conclusion: fitted inhomogeneous Poisson model looks OK

    ######################################################
    # Advanced use of 'expr' argument
    # 
    # set the initial conditions in Metropolis-Hastings algorithm
    # 
    expr <- expression(rmh(fit, start=list(n.start=42), verbose=FALSE))
    \dontrun{qqplot.ppm(fit, 100, expr)}
    \testonly{qqplot.ppm(fit, 4, expr)}

}
\author{Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  and Rolf Turner \email{r.turner@auckland.ac.nz}
  
}
\keyword{spatial}
\keyword{models}
\keyword{hplot}
