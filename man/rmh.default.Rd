\name{rmh.default}
\alias{rmh.default}
\title{Simulate Point Process Models using the Metropolis-Hastings Algorithm.}
\description{
  Generates a random point pattern, simulated from
  a chosen point process model, using the Metropolis-Hastings
  algorithm. 
}

\usage{
   \method{rmh}{default}(model, start=NULL,
   control=default.rmhcontrol(model),
   \dots, 
   verbose=TRUE, snoop=FALSE)
}

\arguments{
  \item{model}{Data specifying the point process model
    that is to be simulated.
  }
  \item{start}{Data determining the initial state of
    the algorithm.
  }
  \item{control}{Data controlling the iterative behaviour
    and termination of the algorithm.
  }
  \item{\dots}{
    Further arguments passed to \code{\link{rmhcontrol}}
    or to trend functions in \code{model}.
  }
  \item{verbose}{Logical flag indicating whether to print progress
    reports.
  }
  \item{snoop}{
    Logical. If \code{TRUE}, activate the visual debugger. 
  }
}

\value{A point pattern (an object of class \code{"ppp"}, see
  \code{\link{ppp.object}}).

  The returned value has an attribute \code{info} containing
  modified versions of the arguments
  \code{model}, \code{start}, and \code{control} which together specify
  the exact simulation procedure. The \code{info} attribute can be
  printed (and is printed automatically by \code{\link{summary.ppp}}).

  The value of \code{\link[base:Random]{.Random.seed}} at the start
  of the simulations is also saved and returned as an attribute
  \code{seed}.

  If the argument \code{track=TRUE} was given (see \code{\link{rmhcontrol}}),
  the transition history of the algorithm
  is saved, and returned as an attribute \code{history}. The transition
  history is a data frame containing a factor \code{proposaltype}
  identifying the proposal type (Birth, Death or Shift) and
  a logical vector \code{accepted} indicating whether the proposal was
  accepted.
  The data frame also has columns \code{numerator}, \code{denominator}
  which give the numerator and denominator of the Hastings ratio for
  the proposal.

  If the argument \code{nsave} was given (see \code{\link{rmhcontrol}}),
  the return value has an attribute \code{saved} which is a list of
  point patterns, containing the intermediate states of the algorithm.
}

\details{
  This function generates simulated realisations from any of a range of
  spatial point processes, using the Metropolis-Hastings algorithm.
  It is the default method for the generic function \code{\link{rmh}}.

  This function executes a Metropolis-Hastings algorithm
  with birth, death and shift proposals as described in
  Geyer and \ifelse{latex}{\out{M\o ller}}{Moller} (1994).

  The argument \code{model} specifies the point process model to be
  simulated. It is either a list, or an object of class
  \code{"rmhmodel"}, with the following components:

  \describe{
    \item{cif}{A character string specifying the choice of
      interpoint interaction for the point process.
    }
    \item{par}{
      Parameter values for the conditional
      intensity function.
    }
    \item{w}{
      (Optional) window in which the pattern is
      to be generated. An object of class \code{"owin"},
      or data acceptable to \code{\link{as.owin}}.
    }
    \item{trend}{
      Data specifying the spatial trend in the model, if it has a trend.
      This may be a function, a pixel image (of class \code{"im"}),
      (or a list of functions or images if the model
      is multitype).
      
      If the trend is a function or functions,
      any auxiliary arguments \code{...} to \code{rmh.default}
      will be passed to these functions, which
      should be of the form \code{function(x, y, ...)}.
    }
    \item{types}{
      List of possible types, for a multitype point process.
    }
  }
  For full details of these parameters, see \code{\link{rmhmodel.default}}.
  
  The argument \code{start} determines the initial state of the
  Metropolis-Hastings algorithm. It is either \code{NULL},
  or an object of class \code{"rmhstart"},
  or a list with the following components:

  \describe{
    \item{n.start}{
      Number of points in the initial point pattern.
      A single integer, or a vector of integers giving the
      numbers of points of each type in a multitype point pattern.
      Incompatible with \code{x.start}.
    }
    \item{x.start}{
      Initial point pattern configuration.
      Incompatible with \code{n.start}.

      \code{x.start} may be a point pattern (an
      object of class \code{"ppp"}), or data which can be coerced
      to this class by \code{\link{as.ppp}},  or an object with
      components \code{x} and \code{y}, or a two-column matrix.
      In the last two cases, the window for the pattern is determined
      by \code{model$w}.
      In the first two cases, if \code{model$w} is also present,
      then the final simulated pattern will be clipped to
      the window \code{model$w}.
    }
  }
  For full details of these parameters, see \code{\link{rmhstart}}.

  The third argument \code{control} controls the simulation
  procedure (including \emph{conditional simulation}),
  iterative behaviour, and termination of the
  Metropolis-Hastings algorithm. It is either \code{NULL}, or
  a list, or an object of class \code{"rmhcontrol"}, with components:
  \describe{
    \item{p}{The probability of proposing a ``shift''
      (as opposed to a birth or death) in the Metropolis-Hastings
      algorithm.
    }
    \item{q}{The conditional probability of proposing a death
      (rather than a birth)
      given that birth/death has been chosen over shift.  
    }
    \item{nrep}{The number of repetitions or iterations
      to be made by the Metropolis-Hastings algorithm.  It should
      be large.
    }
    \item{expand}{
      Either a numerical expansion factor, or
      a window (object of class \code{"owin"}). Indicates that
      the process is to be simulated on a larger domain than the
      original data window \code{w}, then clipped to \code{w}
      when the algorithm has finished.

      The default is to expand the simulation window
      if the model is stationary and non-Poisson
      (i.e. it has no trend and the interaction is not Poisson)
      and not to expand in all other cases. 
      
      If the model has a trend, then in order for expansion to
      be feasible, the trend must be given either as a function,
      or an image whose bounding box is large enough to contain
      the expanded window.
    }
    \item{periodic}{A logical scalar; if \code{periodic} is \code{TRUE}
      we simulate a process on the torus formed by identifying
      opposite edges of a rectangular window.  
    }
    \item{ptypes}{A vector of probabilities (summing to 1) to be used
      in assigning a random type to a new point.
    }
    \item{fixall}{A logical scalar specifying whether to condition on
      the number of points of each type.
    }
    \item{nverb}{An integer specifying how often ``progress reports''
      (which consist simply of the number of repetitions completed)
      should be printed out.  If nverb is left at 0, the default,
      the simulation proceeds silently.
    }
    \item{x.cond}{If this argument is present, then
      \emph{conditional simulation} will be performed, and \code{x.cond}
      specifies the conditioning points and the type of conditioning.
    }
    \item{nsave,nburn}{
      If these values are specified, then
      intermediate states of the simulation algorithm will be saved
      every \code{nsave} iterations, after an initial burn-in period of
      \code{nburn} iterations.
    }
    \item{track}{
      Logical flag indicating whether to save the transition
      history of the simulations.
    }
  }
  For full details of these parameters, see \code{\link{rmhcontrol}}.
  The control parameters can also be given in the \code{\dots} arguments.
}
\section{Conditional Simulation}{
  There are several kinds of conditional simulation.
  \itemize{
    \item
    Simulation \emph{conditional upon the number of points},
    that is, holding the number of points fixed.
    To do this, set \code{control$p} (the probability of a shift) equal to 1.
    The number of points is then determined by the starting state, which
    may be specified either by setting \code{start$n.start} to be a
    scalar, or by setting the initial pattern \code{start$x.start}.
    \item 
    In the case of multitype processes, it is possible to simulate the
    model \emph{conditionally upon the number of points of each type},
    i.e. holding the number of points of each type
    to be fixed. To do this, set \code{control$p} equal to 1
    and \code{control$fixall} to be \code{TRUE}.
    The number of points is then determined by the starting state, which
    may be specified either by setting \code{start$n.start} to be an
    integer vector, or by setting the initial pattern \code{start$x.start}.
    \item
    Simulation 
    \emph{conditional on the configuration observed in a sub-window},
    that is, requiring that, inside a specified sub-window \eqn{V},
    the simulated pattern should agree with a specified point pattern
    \eqn{y}.To do this, set \code{control$x.cond} to equal the
    specified point pattern \eqn{y}, making sure that it is an object of class
    \code{"ppp"} and that the window \code{as.owin(control$x.cond)}
    is the conditioning window \eqn{V}.
    \item
    Simulation \emph{conditional on the presence of specified points},
    that is, requiring that the simulated pattern should include a
    specified set of points. This is simulation from the Palm
    distribution of the point process given a pattern \eqn{y}.
    To do this, set \code{control$x.cond} to be a
    \code{data.frame} containing the coordinates (and marks,
    if appropriate) of the specified points.
  }
  For further information, see \code{\link{rmhcontrol}}.
  
  Note that, when we simulate conditionally on the number of points, or
  conditionally on the number of points of each type,
  no expansion of the window is possible.
}
\section{Visual Debugger}{
  If \code{snoop = TRUE}, an interactive debugger is activated.
  On the current plot device, the debugger displays the current
  state of the Metropolis-Hastings algorithm together with
  the proposed transition to the next state.
  Clicking on this graphical display (using the left mouse button)
  will re-centre the display at the clicked location.
  Surrounding this graphical display is an array of boxes representing
  different actions.
  Clicking on one of the action boxes (using the left mouse button)
  will cause the action to be performed.
  Debugger actions include:
  \itemize{
    \item Zooming in or out
    \item Panning (shifting the field of view) left, right, up or down
    \item Jumping to the next iteration
    \item Skipping 10, 100, 1000, 10000 or 100000 iterations
    \item Jumping to the next Birth proposal (etc)
    \item Changing the fate of the proposal (i.e. changing whether
    the proposal is accepted or rejected)
    \item Dumping the current state and proposal to a file
    \item Printing detailed information at the terminal
    \item Exiting the debugger (so that the simulation
    algorithm continues without further interruption).
  }
  Right-clicking the mouse will also cause the debugger to exit.
}

\references{
   Baddeley, A. and Turner, R. (2000) Practical maximum
   pseudolikelihood for spatial point patterns.
   \emph{Australian and New Zealand Journal of Statistics}
   \bold{42}, 283 -- 322.

   Diggle, P. J. (2003) \emph{Statistical Analysis of Spatial Point
   Patterns} (2nd ed.) Arnold, London.

   Diggle, P.J. and Gratton, R.J. (1984)
   Monte Carlo methods of inference for implicit statistical models.
   \emph{Journal of the Royal Statistical Society, series B}
   \bold{46}, 193 -- 212.

   Diggle, P.J., Gates, D.J., and Stibbard, A. (1987)
   A nonparametric estimator for pairwise-interaction point processes.
   Biometrika \bold{74}, 763 -- 770.

   Geyer, C.J. and \ifelse{latex}{\out{M\o ller}}{Moller}, J. (1994)
   Simulation procedures and likelihood inference for spatial
   point processes.
   \emph{Scandinavian Journal of Statistics} \bold{21}, 359--373.

   Geyer, C.J. (1999)
   Likelihood Inference for Spatial Point
   Processes. Chapter 3 in  O.E. Barndorff-Nielsen, W.S. Kendall and
   M.N.M. Van Lieshout (eds) \emph{Stochastic Geometry: Likelihood and
   Computation}, Chapman and Hall / CRC,  Monographs on Statistics and
   Applied Probability, number 80. Pages 79--140.
}

\section{Warnings}{

There is never a guarantee that the Metropolis-Hastings algorithm
has converged to its limiting distribution.

If \code{start$x.start} is specified then \code{expand} is set equal to 1
and simulation takes place in \code{x.start$window}.  Any specified
value for \code{expand} is simply ignored.

The presence of both a component \code{w} of \code{model} and a
non-null value for \code{x.start$window} makes sense ONLY if \code{w}
is contained in \code{x.start$window}.  

For multitype processes make sure that, even if there is to be no
trend corresponding to a particular type, there is still a component
(a NULL component) for that type, in the list.
}

\seealso{
  \code{\link{rmh}},
  \code{\link{rmh.ppm}},
  \code{\link{rStrauss}},
  \code{\link{ppp}},
  \code{\link{ppm}},
  \code{\link{AreaInter}},
  \code{\link{BadGey}},
  \code{\link{DiggleGatesStibbard}},
  \code{\link{DiggleGratton}},
  \code{\link{Fiksel}},
  \code{\link{Geyer}},
  \code{\link{Hardcore}},
  \code{\link{LennardJones}},
  \code{\link{MultiHard}},
  \code{\link{MultiStrauss}},
  \code{\link{MultiStraussHard}},
  \code{\link{PairPiece}},
  \code{\link{Poisson}},
  \code{\link{Softcore}},
  \code{\link{Strauss}},
  \code{\link{StraussHard}},
  \code{\link{Triplets}}
}
\section{Other models}{
  In theory, any finite point process model can be simulated using
  the Metropolis-Hastings algorithm, provided the conditional
  intensity is uniformly bounded.

  In practice, the list of point process models that can be simulated using
  \code{rmh.default} is limited to those that have been implemented
  in the package's internal C code. More options will be added in the future.

  Note that the \code{lookup} conditional intensity function
  permits the simulation (in theory, to any desired degree
  of approximation) of any pairwise interaction process for
  which the interaction depends only on the distance between
  the pair of points.
}
\section{Reproducible simulations}{
  If the user wants the simulation to be exactly reproducible
  (e.g. for a figure in a journal article, where it is useful to
  have the figure consistent from draft to draft) then the state of
  the random number generator should be set before calling
  \code{rmh.default}. This can be done either by calling
  \code{\link[base:Random]{set.seed}} or by assigning a value to
  \code{\link[base:Random]{.Random.seed}}. In the examples below, we use
  \code{\link[base:Random]{set.seed}}.  

  If a simulation has been performed and the user now wants to 
  repeat it exactly, the random seed should be extracted from
  the simulated point pattern \code{X} by \code{seed <- attr(x, "seed")},
  then assigned to the system random nunber state by
  \code{.Random.seed <- seed} before calling \code{rmh.default}.
}
\examples{
   if(interactive()) {
     nr   <- 1e5
     nv  <- 5000
     ns <- 200
   } else {
     nr  <- 10
     nv <- 5
     ns <- 20
     oldopt <- spatstat.options()
     spatstat.options(expand=1.1)
   }
   set.seed(961018)
   
   # Strauss process.
   mod01 <- list(cif="strauss",par=list(beta=2,gamma=0.2,r=0.7),
                 w=c(0,10,0,10))
   X1.strauss <- rmh(model=mod01,start=list(n.start=ns),
                     control=list(nrep=nr,nverb=nv))

   if(interactive()) plot(X1.strauss)
   
   # Strauss process, conditioning on n = 42:
   X2.strauss <- rmh(model=mod01,start=list(n.start=42),
                     control=list(p=1,nrep=nr,nverb=nv))

   # Tracking algorithm progress:
   X <- rmh(model=mod01,start=list(n.start=ns),
            control=list(nrep=nr, nsave=nr/5, nburn=nr/2, track=TRUE))
   History <- attr(X, "history")
   Saved <- attr(X, "saved")
   head(History)
   plot(Saved)

   # Hard core process:
   mod02 <- list(cif="hardcore",par=list(beta=2,hc=0.7),w=c(0,10,0,10))
   X3.hardcore <- rmh(model=mod02,start=list(n.start=ns),
                     control=list(nrep=nr,nverb=nv))
   
   if(interactive()) plot(X3.hardcore)

   # Strauss process equal to pure hardcore:
   mod02s <- list(cif="strauss",par=list(beta=2,gamma=0,r=0.7),w=c(0,10,0,10))
   X3.strauss <- rmh(model=mod02s,start=list(n.start=ns),
                     control=list(nrep=nr,nverb=nv))
   
   # Strauss process in a polygonal window.
   x     <- c(0.55,0.68,0.75,0.58,0.39,0.37,0.19,0.26,0.42)
   y     <- c(0.20,0.27,0.68,0.99,0.80,0.61,0.45,0.28,0.33)
   mod03 <- list(cif="strauss",par=list(beta=2000,gamma=0.6,r=0.07),
                w=owin(poly=list(x=x,y=y)))
   X4.strauss <- rmh(model=mod03,start=list(n.start=ns),
                     control=list(nrep=nr,nverb=nv))
   if(interactive()) plot(X4.strauss)
   
   # Strauss process in a polygonal window, conditioning on n = 80.
   X5.strauss <- rmh(model=mod03,start=list(n.start=ns),
                     control=list(p=1,nrep=nr,nverb=nv))
   
   # Strauss process, starting off from X4.strauss, but with the
   # polygonal window replace by a rectangular one.  At the end,
   # the generated pattern is clipped to the original polygonal window.
   xxx <- X4.strauss
   xxx$window <- as.owin(c(0,1,0,1))
   X6.strauss <- rmh(model=mod03,start=list(x.start=xxx),
                     control=list(nrep=nr,nverb=nv))
   
   # Strauss with hardcore:
   mod04 <- list(cif="straush",par=list(beta=2,gamma=0.2,r=0.7,hc=0.3),
                w=c(0,10,0,10))
   X1.straush <- rmh(model=mod04,start=list(n.start=ns),
                     control=list(nrep=nr,nverb=nv))
   
   # Another Strauss with hardcore (with a perhaps surprising result):
   mod05 <- list(cif="straush",par=list(beta=80,gamma=0.36,r=45,hc=2.5),
                w=c(0,250,0,250))
   X2.straush <- rmh(model=mod05,start=list(n.start=ns),
                     control=list(nrep=nr,nverb=nv))
   
   # Pure hardcore (identical to X3.strauss).
   mod06 <- list(cif="straush",par=list(beta=2,gamma=1,r=1,hc=0.7),
                w=c(0,10,0,10))
   X3.straush <- rmh(model=mod06,start=list(n.start=ns),
                     control=list(nrep=nr,nverb=nv))
   
   # Soft core:
   w    <- c(0,10,0,10)
   mod07 <- list(cif="sftcr",par=list(beta=0.8,sigma=0.1,kappa=0.5),
                w=c(0,10,0,10))
   X.sftcr <- rmh(model=mod07,start=list(n.start=ns),
                  control=list(nrep=nr,nverb=nv))
   if(interactive()) plot(X.sftcr)

   # Area-interaction process:
   mod42 <- rmhmodel(cif="areaint",par=list(beta=2,eta=1.6,r=0.7),
                 w=c(0,10,0,10))
   X.area <- rmh(model=mod42,start=list(n.start=ns),
                  control=list(nrep=nr,nverb=nv))
   if(interactive()) plot(X.area)

   # Triplets process
   modtrip <- list(cif="triplets",par=list(beta=2,gamma=0.2,r=0.7),
                   w=c(0,10,0,10))
   X.triplets <- rmh(model=modtrip,
                     start=list(n.start=ns),
                     control=list(nrep=nr,nverb=nv))
   if(interactive()) plot(X.triplets)
   
   # Multitype Strauss:
   beta <- c(0.027,0.008)
   gmma <- matrix(c(0.43,0.98,0.98,0.36),2,2)
   r    <- matrix(c(45,45,45,45),2,2)
   mod08 <- list(cif="straussm",par=list(beta=beta,gamma=gmma,radii=r),
                w=c(0,250,0,250))
   X1.straussm <- rmh(model=mod08,start=list(n.start=ns),
                      control=list(ptypes=c(0.75,0.25),nrep=nr,nverb=nv))
   if(interactive()) plot(X1.straussm)
   
   # Multitype Strauss conditioning upon the total number
   # of points being 80:
   X2.straussm <- rmh(model=mod08,start=list(n.start=ns),
                      control=list(p=1,ptypes=c(0.75,0.25),nrep=nr,
                                   nverb=nv))
   
   # Conditioning upon the number of points of type 1 being 60
   # and the number of points of type 2 being 20:
   X3.straussm <- rmh(model=mod08,start=list(n.start=c(60,20)),
                      control=list(fixall=TRUE,p=1,ptypes=c(0.75,0.25),
                                   nrep=nr,nverb=nv))
   
   # Multitype Strauss hardcore:
   rhc  <- matrix(c(9.1,5.0,5.0,2.5),2,2)
   mod09 <- list(cif="straushm",par=list(beta=beta,gamma=gmma,
                iradii=r,hradii=rhc),w=c(0,250,0,250))
   X.straushm <- rmh(model=mod09,start=list(n.start=ns),
                     control=list(ptypes=c(0.75,0.25),nrep=nr,nverb=nv))
   
   # Multitype Strauss hardcore with trends for each type:
   beta  <- c(0.27,0.08)
   tr3   <- function(x,y){x <- x/250; y <- y/250;
   			   exp((6*x + 5*y - 18*x^2 + 12*x*y - 9*y^2)/6)
                         }
                         # log quadratic trend
   tr4   <- function(x,y){x <- x/250; y <- y/250;
                         exp(-0.6*x+0.5*y)}
                        # log linear trend
   mod10 <- list(cif="straushm",par=list(beta=beta,gamma=gmma,
                 iradii=r,hradii=rhc),w=c(0,250,0,250),
                 trend=list(tr3,tr4))
   X1.straushm.trend <- rmh(model=mod10,start=list(n.start=ns),
                            control=list(ptypes=c(0.75,0.25),
                            nrep=nr,nverb=nv))
   if(interactive()) plot(X1.straushm.trend)
   
   # Multitype Strauss hardcore with trends for each type, given as images:
   bigwin <- square(250)
   i1 <- as.im(tr3, bigwin)
   i2 <- as.im(tr4, bigwin)   
   mod11 <- list(cif="straushm",par=list(beta=beta,gamma=gmma,
                 iradii=r,hradii=rhc),w=bigwin,
                 trend=list(i1,i2))
   X2.straushm.trend <- rmh(model=mod11,start=list(n.start=ns),
                            control=list(ptypes=c(0.75,0.25),expand=1,
                            nrep=nr,nverb=nv))
   
   # Diggle, Gates, and Stibbard:
   mod12 <- list(cif="dgs",par=list(beta=3600,rho=0.08),w=c(0,1,0,1))
   X.dgs <- rmh(model=mod12,start=list(n.start=ns),
                control=list(nrep=nr,nverb=nv))
   if(interactive()) plot(X.dgs)
   
   # Diggle-Gratton:
   mod13 <- list(cif="diggra",
                 par=list(beta=1800,kappa=3,delta=0.02,rho=0.04),
                 w=square(1))
   X.diggra <- rmh(model=mod13,start=list(n.start=ns),
                   control=list(nrep=nr,nverb=nv))
   if(interactive()) plot(X.diggra)
   
   # Fiksel:
   modFik <- list(cif="fiksel",
                 par=list(beta=180,r=0.15,hc=0.07,kappa=2,a= -1.0),
                 w=square(1))
   X.fiksel <- rmh(model=modFik,start=list(n.start=ns),
                   control=list(nrep=nr,nverb=nv))
   if(interactive()) plot(X.fiksel)
   
   # Geyer:
   mod14 <- list(cif="geyer",par=list(beta=1.25,gamma=1.6,r=0.2,sat=4.5),
                 w=c(0,10,0,10))
   X1.geyer <- rmh(model=mod14,start=list(n.start=ns),
                   control=list(nrep=nr,nverb=nv))
   if(interactive()) plot(X1.geyer)
   
   # Geyer; same as a Strauss process with parameters
   # (beta=2.25,gamma=0.16,r=0.7):
   
   mod15 <- list(cif="geyer",par=list(beta=2.25,gamma=0.4,r=0.7,sat=10000),
                 w=c(0,10,0,10))
   X2.geyer <- rmh(model=mod15,start=list(n.start=ns),
                   control=list(nrep=nr,nverb=nv))
   
   mod16 <- list(cif="geyer",par=list(beta=8.1,gamma=2.2,r=0.08,sat=3))
   data(redwood)
   X3.geyer <- rmh(model=mod16,start=list(x.start=redwood),
                   control=list(periodic=TRUE,nrep=nr,nverb=nv))
   
   # Geyer, starting from the redwood data set, simulating
   # on a torus, and conditioning on n:
   X4.geyer <- rmh(model=mod16,start=list(x.start=redwood),
                   control=list(p=1,periodic=TRUE,nrep=nr,nverb=nv))

   # Lookup (interaction function h_2 from page 76, Diggle (2003)):
      r <- seq(from=0,to=0.2,length=101)[-1] # Drop 0.
      h <- 20*(r-0.05)
      h[r<0.05] <- 0
      h[r>0.10] <- 1
      mod17 <- list(cif="lookup",par=list(beta=4000,h=h,r=r),w=c(0,1,0,1))
      X.lookup <- rmh(model=mod17,start=list(n.start=ns),
                      control=list(nrep=nr,nverb=nv))
      if(interactive()) plot(X.lookup)
                   
   # Strauss with trend
   tr <- function(x,y){x <- x/250; y <- y/250;
   			   exp((6*x + 5*y - 18*x^2 + 12*x*y - 9*y^2)/6)
                         }
   beta <- 0.3
   gmma <- 0.5
   r    <- 45
   modStr <- list(cif="strauss",par=list(beta=beta,gamma=gmma,r=r),
                 w=square(250), trend=tr)
   X1.strauss.trend <- rmh(model=modStr,start=list(n.start=ns),
                           control=list(nrep=nr,nverb=nv))
   # Baddeley-Geyer
   r <- seq(0,0.2,length=8)[-1]
   gmma <- c(0.5,0.6,0.7,0.8,0.7,0.6,0.5)
   mod18 <- list(cif="badgey",par=list(beta=4000, gamma=gmma,r=r,sat=5),
                 w=square(1))
   X1.badgey <- rmh(model=mod18,start=list(n.start=ns),
                    control=list(nrep=nr,nverb=nv))
   mod19 <- list(cif="badgey",
                 par=list(beta=4000, gamma=gmma,r=r,sat=1e4),
                 w=square(1))
   set.seed(1329)
   X2.badgey <- rmh(model=mod18,start=list(n.start=ns),
                    control=list(nrep=nr,nverb=nv))

   # Check:
   h <- ((prod(gmma)/cumprod(c(1,gmma)))[-8])^2
   hs <- stepfun(r,c(h,1))
   mod20 <- list(cif="lookup",par=list(beta=4000,h=hs),w=square(1))
   set.seed(1329)
   X.check <- rmh(model=mod20,start=list(n.start=ns),
                      control=list(nrep=nr,nverb=nv))
   # X2.badgey and X.check will be identical.

   mod21 <- list(cif="badgey",par=list(beta=300,gamma=c(1,0.4,1),
                 r=c(0.035,0.07,0.14),sat=5), w=square(1))
   X3.badgey <- rmh(model=mod21,start=list(n.start=ns),
                    control=list(nrep=nr,nverb=nv))
   # Same result as Geyer model with beta=300, gamma=0.4, r=0.07,
   # sat = 5 (if seeds and control parameters are the same)

   # Or more simply:
   mod22 <- list(cif="badgey",
                 par=list(beta=300,gamma=0.4,r=0.07, sat=5),
                 w=square(1))
   X4.badgey <- rmh(model=mod22,start=list(n.start=ns),
                    control=list(nrep=nr,nverb=nv))
   # Same again --- i.e. the BadGey model includes the Geyer model.


   # Illustrating scalability.
   \dontrun{
    M1 <- rmhmodel(cif="strauss",par=list(beta=60,gamma=0.5,r=0.04),w=owin())
    set.seed(496)
    X1 <- rmh(model=M1,start=list(n.start=300))
    M2 <- rmhmodel(cif="strauss",par=list(beta=0.6,gamma=0.5,r=0.4),
              w=owin(c(0,10),c(0,10)))
    set.seed(496)
    X2  <- rmh(model=M2,start=list(n.start=300))
    chk <- affine(X1,mat=diag(c(10,10)))
    all.equal(chk,X2,check.attributes=FALSE)
    # Under the default spatstat options the foregoing all.equal()
    # will yield TRUE.  Setting spatstat.options(scalable=FALSE) and
    # re-running the code will reveal differences between X1 and X2.
   }

   if(!interactive()) spatstat.options(oldopt)
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{datagen}
