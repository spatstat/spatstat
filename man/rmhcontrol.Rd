\name{rmhcontrol}
\alias{rmhcontrol}
\alias{rmhcontrol.default}
\title{Set Control Parameters for Metropolis-Hastings Algorithm.}
\description{
  Sets up a list of parameters controlling the iterative behaviour
  of the Metropolis-Hastings algorithm. 
}
\usage{
   rmhcontrol(\dots)

   \method{rmhcontrol}{default}(\dots, p=0.9, q=0.5, nrep=5e5,
                      expand=NULL, periodic=NULL, ptypes=NULL,
                      x.cond=NULL, fixall=FALSE, nverb=0,
                      nsave=NULL, nburn=nsave, track=FALSE)
}
\arguments{
  \item{\dots}{Arguments passed to methods.}
  \item{p}{Probability of proposing a shift (as against a birth/death).}
  \item{q}{Conditional probability of proposing a death given that a
    birth or death will be proposed.}
  \item{nrep}{Total number of steps (proposals) of Metropolis-Hastings
    algorithm that should be run.}
  \item{expand}{
    Simulation window or expansion rule.
    Either a window (object of class \code{"owin"})
    or a numerical expansion factor, specifying that
    simulations are to be performed in a domain other than the
    original data window, then clipped to the original data window.
    This argument is passed to \code{\link{rmhexpand}}.
    A numerical expansion factor can be in several formats:
    see \code{\link{rmhexpand}}.
  }
  \item{periodic}{
    Logical value (or \code{NULL}) indicating whether to simulate
    ``periodically'', i.e. identifying opposite edges of the rectangular
    simulation window. A \code{NULL} value means ``undecided.''
  }
  \item{ptypes}{For multitype point processes, the distribution of
    the mark attached to a new random point (when a birth is
    proposed)}
  \item{x.cond}{Conditioning points for conditional simulation.}
  \item{fixall}{(Logical) for multitype point processes,
    whether to fix the number of points of each type.}
  \item{nverb}{Progress reports will be printed every \code{nverb}
    iterations}
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

\value{
  An object of class \code{"rmhcontrol"}, which is essentially
  a list of parameter values for the algorithm.
  
  There is a \code{print} method for this class, which prints
  a sensible description of the parameters chosen.
}
\details{
  The Metropolis-Hastings algorithm, implemented as \code{\link{rmh}},
  generates simulated realisations of point process models.
  The function \code{rmhcontrol}
  sets up a list of parameters which control the 
  iterative behaviour
  and termination of the Metropolis-Hastings algorithm, for use in a
  subsequent call to \code{\link{rmh}}. It also checks that the
  parameters are valid.

  (A separate function \code{\link{rmhstart}}
  determines the initial state of the algorithm,
  and \code{\link{rmhmodel}} determines the model to be simulated.)

  The parameters are as follows:
  
  \describe{
    \item{p}{The probability of proposing a ``shift''
      (as opposed to a birth or death) in the Metropolis-Hastings
      algorithm.

      If \eqn{p = 1} then the algorithm only alters existing points,
      so the number of points never changes, i.e. we are
      simulating conditionally upon the number of points.
      The number of points is determined by the initial state
      (specified by \code{\link{rmhstart}}).

      If \eqn{p=1} and \code{fixall=TRUE} and the model
      is a multitype point process model, then the algorithm
      only shifts the locations of existing points and does not
      alter their marks (types). 
      This is equivalent to simulating conditionally
      upon the number of points of each type.
      These numbers are again specified by the initial state.

      If \eqn{p = 1} then no expansion of the simulation window
      is allowed (see \code{expand} below).

      The default value of \code{p} can be changed by setting
      the parameter \code{rmh.p} in \code{\link{spatstat.options}}.
    }
    \item{q}{The conditional probability of proposing a death
      (rather than a birth) given that a shift is not proposed.
      This is of course ignored if \code{p} is equal to 1.

      The default value of \code{q} can be changed by setting
      the parameter \code{rmh.q} in \code{\link{spatstat.options}}.
    }
    \item{nrep}{The number of repetitions or iterations 
      to be made by the Metropolis-Hastings algorithm.  It should
      be large.

      The default value of \code{nrep} can be changed by setting
      the parameter \code{rmh.nrep} in \code{\link{spatstat.options}}.
    }
    \item{expand}{
      Either a number or a window (object of class \code{"owin"}).
      Indicates that the process is to be simulated on a 
      domain other than the original data window \code{w},
      then clipped to \code{w} when the algorithm has finished.
      This would often be done in order to approximate the
      simulation of a stationary process (Geyer, 1999)
      or more generally a process existing in the
      whole plane, rather than just in the window \code{w}.

      If \code{expand} is a window object, it is taken as the
      larger domain in which simulation is performed.

      If \code{expand} is numeric, it is interpreted
      as an expansion factor or expansion distance
      for determining the simulation domain from the data window.
      It should be a \emph{named} scalar, such as
      \code{expand=c(area=2)}, \code{expand=c(distance=0.1)},
      \code{expand=c(length=1.2)}.  See \code{\link{rmhexpand}()} for
      more details. If the name is omitted, it defaults to \code{area}.
     
      Expansion is not permitted if the number of points has been
      fixed by setting \code{p = 1} or if the
      starting configuration has been specified via the
      argument \code{x.start} in \code{\link{rmhstart}}.

      If \code{expand} is \code{NULL}, this is interpreted to mean
      \dQuote{not yet decided}. An expansion rule will be determined
      at a later stage, using appropriate defaults.
      See \code{\link{rmhexpand}}.
    }
    \item{periodic}{A logical value (or \code{NULL})
      determining whether to simulate \dQuote{periodically}.
      If \code{periodic} is \code{TRUE}, and if the simulation window
      is a rectangle, then the simulation algorithm effectively
      identifies opposite edges of the rectangle. Points
      near the right-hand edge of the rectangle are deemed to be close
      to points near the left-hand edge. Periodic simulation usually
      gives a better approximation to a stationary point process.
      For periodic simulation, the simulation window must be a rectangle.
      (The simulation window is determined by \code{expand} as described
      above.)

      The value \code{NULL} means \sQuote{undecided}.
      The decision is postponed until \code{\link{rmh}} is called.
      Depending on the point process model to be simulated,
      \code{rmh} will then set \code{periodic=TRUE} if the simulation window
      is expanded \emph{and} the expanded simulation window is rectangular;
      otherwise \code{periodic=FALSE}.

      Note that \code{periodic=TRUE} is only permitted when the
      simulation window (i.e. the expanded window) is rectangular.       
    }
    \item{ptypes}{A vector of probabilities (summing to 1) to be used
      in assigning a random type to a new point.  Defaults to a vector
      each of whose entries is \eqn{1/nt} where \eqn{nt} is the number
      of types for the process.  Convergence of the simulation
      algorithm should be improved if \code{ptypes} is close to the
      relative frequencies of the types which will result from the
      simulation.
    }
    \item{x.cond}{
      If this argument is given,
      then \emph{conditional simulation} will be performed,
      and \code{x.cond} specifies the location of the
      fixed points as well as the type of conditioning.
      It should be either a point pattern
      (object of class \code{"ppp"}) or a \code{list(x,y)}
      or a \code{data.frame}.
      See the section on Conditional Simulation.
    }
    \item{fixall}{A logical scalar specifying whether to condition on
      the number of points of each type.  Meaningful only if a marked
      process is being simulated, and if \eqn{p = 1}.  A warning message
      is given if \code{fixall} is set equal to \code{TRUE} when it is
      not meaningful.
    }
    \item{nverb}{An integer specifying how often ``progress reports''
      (which consist simply of the number of repetitions completed)
      should be printed out.  If nverb is left at 0, the default,
      the simulation proceeds silently.
    }
    \item{nsave,nburn}{If these integers are given, then the
      current state of the simulation algorithm (i.e. the current
      random point pattern) will be saved every \code{nsave} iterations,
      starting from iteration \code{nburn}. 
    }
    \item{track}{
      Logical flag indicating whether to save the transition
      history of the simulations (i.e. information specifying
      what type of proposal was made, and whether it was accepted
      or rejected, for each iteration).
    }
  }
}
\section{Conditional Simulation}{
  For a Gibbs point process \eqn{X},
  the Metropolis-Hastings algorithm easily accommodates several
  kinds of conditional simulation:
  \describe{
    \item{conditioning on the total number of points:}{
      We fix the total number of points \eqn{N(X)} to be equal to
      \eqn{n}. We simulate from the conditional distribution of
      \eqn{X} given \eqn{N(X) = n}. 
    }
    \item{conditioning on the number of points of each type:}{
      In a multitype point process, where \eqn{Y_j}{Y[[j]]} denotes the
      process of points of type \eqn{j}, we fix the number
      \eqn{N(Y_j)}{N(Y[[j]])} of points of type \eqn{j} to be equal to
      \eqn{n_j}{n[j]}, for \eqn{j=1,2,\ldots,m}{j=1,2,...,m}.
      We simulate from the conditional distribution of \eqn{X}
      given \eqn{N(Y_j)=n_j}{N(Y[[j]]) = n[j]} for
      \eqn{j=1,2,\ldots,m}{j=1,2,...,m}. 
    }
    \item{conditioning on the realisation in a subwindow:}{
      We require that the point process \eqn{X} should,
      within a specified sub-window \eqn{V},
      coincide with a specified point pattern \eqn{y}.
      We simulate from the conditional distribution of \eqn{X}
      given \eqn{X \cap V = y}{(X intersect V) = y}.
    }
    \item{Palm conditioning:}{
      We require that the point process \eqn{X} include
      a specified list of points \eqn{y}. We simulate from
      the point process with probability density
      \eqn{g(x) = c f(x \cup y)}{g(x) = c * f(x union y)}
      where \eqn{f} is the probability density of the original
      process \eqn{X}, and \eqn{c} is a normalising constant.
    }
  }
  To achieve each of these types of conditioning we
  do as follows:
  \describe{
    \item{conditioning on the total number of points:}{
      Set \code{p=1}.
      The number of points is determined by the initial state
      of the simulation: see \code{\link{rmhstart}}.
    }
    \item{conditioning on the number of points of each type:}{
      Set \code{p=1} and \code{fixall=TRUE}.
      The number of points of each type is determined by the initial state
      of the simulation: see \code{\link{rmhstart}}.
    }
    \item{conditioning on the realisation in a subwindow:}{
      Set \code{x.cond} to be a point pattern (object of
      class \code{"ppp"}). Its window \code{V=x.cond$window} becomes the
      conditioning subwindow \eqn{V}.
    }
    \item{Palm conditioning:}{
      Set \code{x.cond} to be a \code{list(x,y)} or \code{data.frame}
      with two columns containing the coordinates of the points, or a 
      \code{list(x,y,marks)} or \code{data.frame} with three columns
      containing the coordinates and marks of the points.
    }
  }
  The arguments \code{x.cond}, \code{p} and \code{fixall} can be
  combined.   
}
\references{
   Geyer, C.J. (1999)
   Likelihood Inference for Spatial Point
   Processes. Chapter 3 in  O.E. Barndorff-Nielsen, W.S. Kendall and
   M.N.M. Van Lieshout (eds) \emph{Stochastic Geometry: Likelihood and
   Computation}, Chapman and Hall / CRC,  Monographs on Statistics and
   Applied Probability, number 80. Pages 79--140.
}

\seealso{
  \code{\link{rmh}},
  \code{\link{rmhmodel}},
  \code{\link{rmhstart}},
  \code{\link{rmhexpand}},
  \code{\link{spatstat.options}}
}
\examples{
   # parameters given as named arguments
   c1 <- rmhcontrol(p=0.3,periodic=TRUE,nrep=1e6,nverb=1e5)

   # parameters given as a list
   liz <- list(p=0.9, nrep=1e4)
   c2 <- rmhcontrol(liz)

   # parameters given in rmhcontrol object
   c3 <- rmhcontrol(c1)
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{datagen}
