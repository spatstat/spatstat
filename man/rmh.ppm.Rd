\name{rmh.ppm}
\alias{rmh.ppm}
\title{Simulate from a Fitted Point Process Model}
\description{
  Given a point process model fitted to data, 
  generate a random simulation of the model, 
  using the Metropolis-Hastings algorithm.
}

\usage{
  \method{rmh}{ppm}(model, start=NULL,
                    control=default.rmhcontrol(model),
                    \dots,
                    project=TRUE, verbose=TRUE)
}

\arguments{
  \item{model}{A fitted point process model (object of class
    \code{"ppm"}, see \code{\link{ppm.object}}) which it is desired
    to simulate.  This fitted model is usually the result of a call
    to \code{\link{ppm}}.  See \bold{Details} below.
  }
  \item{start}{Data determining the initial state
    of the Metropolis-Hastings algorithm.  See
    \code{\link{rmhstart}} for description of these arguments.
    Defaults to \code{list(x.start=data.ppm(model))}
  }
  \item{control}{Data controlling the iterative behaviour of
    the Metropolis-Hastings algorithm.  See \code{\link{rmhcontrol}}
    for description of these arguments.
  }
  \item{\dots}{
    Further arguments passed to \code{\link{rmhcontrol}},
    or to \code{\link{rmh.default}}, or to covariate functions in the model.
  }
  \item{project}{
    Logical flag indicating what to do if the fitted model is
    invalid (in the sense that the values of the fitted coefficients do not
    specify a valid point process).
    If \code{project=TRUE} the closest valid model will be simulated;
    if \code{project=FALSE} an error will occur.
  }
  \item{verbose}{
    Logical flag indicating whether to print progress reports.
  }
}

\value{A point pattern (an object of class \code{"ppp"}; see
  \code{\link{ppp.object}}).
}

\details{
  This function generates simulated realisations from a point
  process model that has been fitted to point pattern data.  It is
  a method for the generic function \code{\link{rmh}} for the
  class \code{"ppm"} of fitted point process models.  To simulate
  other kinds of point process models, see \code{\link{rmh}}
  or \code{\link{rmh.default}}.

  The argument \code{model} describes the fitted model.  It must be
  an object of class \code{"ppm"} (see \code{\link{ppm.object}}),
  and will typically be the result of a call to the point process
  model fitting function \code{\link{ppm}}.

  The current implementation enables simulation from any fitted model
  involving the interactions
  \code{\link{AreaInter}},
  \code{\link{DiggleGratton}},
  \code{\link{DiggleGatesStibbard}},
  \code{\link{Geyer}},
  \code{\link{Hardcore}},
  \code{\link{MultiStrauss}},
  \code{\link{MultiStraussHard}},
  \code{\link{PairPiece}},
  \code{\link{Poisson}},
  \code{\link{Strauss}},
  \code{\link{StraussHard}}
  and \code{\link{Softcore}},
  including nonstationary models. See the examples.
  
  It is also possible to simulate \emph{hybrids} of several such models.
  See \code{\link{Hybrid}} and the examples.
  
  It is possible that the fitted coefficients of a point process model
  may be ``illegal'', i.e. that there may not exist a
  mathematically well-defined point process with the given parameter
  values. For example, a Strauss process with interaction
  parameter \eqn{\gamma > 1}{gamma > 1} does not exist,
  but the model-fitting procedure used in \code{\link{ppm}} will sometimes
  produce values of \eqn{\gamma}{gamma} greater than 1.
  In such cases, if \code{project=FALSE} then an error will occur,
  while if \code{project=TRUE} then \code{rmh.ppm} will find
  the nearest legal model and simulate
  this model instead. (The nearest legal model is obtained by
  projecting the vector of coefficients onto the set of
  valid coefficient vectors. The result is usually the Poisson process
  with the same fitted intensity.)
  
  The arguments \code{start} and \code{control} are lists of
  parameters determining the initial state and the iterative
  behaviour, respectively, of the Metropolis-Hastings algorithm.

  The argument \code{start} is passed directly to \code{\link{rmhstart}}.
  See \code{\link{rmhstart}} for details of the parameters of the
  initial state, and their default values.

  The argument \code{control} is first passed to
  \code{\link{rmhcontrol}}. Then if any additional arguments \code{\dots}
  are given, \code{\link{update.rmhcontrol}} is called to update the
  parameter values. See \code{\link{rmhcontrol}} for details of
  the iterative behaviour parameters, and \code{\link{default.rmhcontrol}}
  for their default values.

  Note that if you specify expansion of the simulation window
  using the parameter \code{expand} (so that the
  model will be simulated on a window larger than the original data
  window) then the model must be capable of extrapolation to this
  larger window. This is usually not possible for models which
  depend on external covariates, because the domain of a covariate image
  is usually the same as the domain of the fitted model.
  
  After extracting the relevant information from the fitted model
  object \code{model}, \code{rmh.ppm} invokes the default
  \code{rmh} algorithm \code{\link{rmh.default}}, unless the model
  is Poisson. If the model is Poisson then the Metropolis-Hastings
  algorithm is not needed, and the model is simulated directly, using
  one of \code{\link{rpoispp}}, \code{\link{rmpoispp}},
  \code{\link{rpoint}} or \code{\link{rmpoint}}.

  See \code{\link{rmh.default}} for further information about the
  implementation, or about the Metropolis-Hastings algorithm.
}

\section{Warnings}{
  See Warnings in \code{\link{rmh.default}}.
}

\seealso{
  \code{\link{simulate.ppm}},
  \code{\link{rmh}},
  \code{\link{rmhmodel}},
  \code{\link{rmhcontrol}},
  \code{\link{default.rmhcontrol}},
  \code{\link{update.rmhcontrol}},
  \code{\link{rmhstart}},
  \code{\link{rmh.default}},
  \code{\link{ppp.object}},
  \code{\link{ppm}},

  Interactions:
  \code{\link{AreaInter}},
  \code{\link{DiggleGratton}},
  \code{\link{DiggleGatesStibbard}},
  \code{\link{Geyer}},
  \code{\link{Hardcore}},
  \code{\link{Hybrid}},
  \code{\link{MultiStrauss}},
  \code{\link{MultiStraussHard}},
  \code{\link{PairPiece}},
  \code{\link{Poisson}},
  \code{\link{Strauss}},
  \code{\link{StraussHard}},
  \code{\link{Softcore}}
}

\examples{
   live <- interactive()
   op <- spatstat.options()
   spatstat.options(rmh.nrep=1e5)
   Nrep <- 1e5

   X <- swedishpines
   if(live) plot(X, main="Swedish Pines data")

   # Poisson process
   fit <- ppm(X, ~1, Poisson())
   Xsim <- rmh(fit)
   if(live) plot(Xsim, main="simulation from fitted Poisson model")

   # Strauss process   
   fit <- ppm(X, ~1, Strauss(r=7))
   Xsim <- rmh(fit)
   if(live) plot(Xsim, main="simulation from fitted Strauss model")

   \dontrun{
     # Strauss process simulated on a larger window
     # then clipped to original window
     Xsim <- rmh(fit, control=list(nrep=Nrep, expand=1.1, periodic=TRUE))
     Xsim <- rmh(fit, nrep=Nrep, expand=2, periodic=TRUE)
   }

   \dontrun{
     X <- rSSI(0.05, 100)
     # piecewise-constant pairwise interaction function
     fit <- ppm(X, ~1, PairPiece(seq(0.02, 0.1, by=0.01)))
     Xsim <- rmh(fit)
   }

    # marked point pattern
    Y <- amacrine

   \dontrun{
     # marked Poisson models
     fit <- ppm(Y)
     fit <- ppm(Y,~marks)
     fit <- ppm(Y,~polynom(x,2))
     fit <- ppm(Y,~marks+polynom(x,2))
     fit <- ppm(Y,~marks*polynom(x,y,2))
     Ysim <- rmh(fit)
   }

   # multitype Strauss models
   MS <- MultiStrauss(types = levels(Y$marks),
                      radii=matrix(0.07, ncol=2, nrow=2))
   \dontrun{
    fit <- ppm(Y, ~marks, MS)
    Ysim <- rmh(fit)
   }

   fit <- ppm(Y,~marks*polynom(x,y,2), MS)
   Ysim <- rmh(fit)
   if(live) plot(Ysim, main="simulation from fitted inhomogeneous Multitype Strauss")

   spatstat.options(op)

  \dontrun{
    # Hybrid model
    fit <- ppm(redwood, ~1, Hybrid(A=Strauss(0.02), B=Geyer(0.1, 2)))
    Y <- rmh(fit)
  }
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{models}
\keyword{datagen}
