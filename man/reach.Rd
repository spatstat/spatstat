\name{reach}
\alias{reach}
\alias{reach.ppm}
\alias{reach.interact}
\alias{reach.rmhmodel}
\title{Interaction Distance of a Point Process}
\description{
  Computes the interaction distance of a point process.
}
\usage{
  reach(x, \dots)
  \method{reach}{ppm}(x, \dots, epsilon=0)
  \method{reach}{interact}(x, \dots)
  \method{reach}{rmhmodel}(x, \dots)
}
\arguments{
  \item{x}{Either a fitted point process model (object of class
    \code{"ppm"}), an interpoint interaction (object of class
    \code{"interact"}) or a point process model for simulation
    (object of class \code{"rmhmodel"}).
  }
  \item{epsilon}{
    Numerical threshold below which interaction is treated as zero.
    See details.
  }
  \item{\dots}{
    Other arguments are ignored.
  }
}
\value{
  The interaction distance, or \code{NA} if this cannot be
  computed from the information given.
}
\details{
  The `interaction distance' or `interaction range' of a point process
  model is the smallest distance \eqn{D} such that any two points in the
  process which are separated by a distance greater than \eqn{D} do not
  interact with each other.

  For example, the interaction range of a Strauss process
  (see \code{\link{Strauss}})
  with parameters \eqn{\beta,\gamma,r}{beta,gamma,r} is equal to
  \eqn{r}, unless \eqn{\gamma=1}{gamma=1} in which case the model is
  Poisson and the interaction
  range is \eqn{0}.
  The interaction range of a Poisson process is zero.
  The interaction range of the Ord threshold process
  (see \code{\link{OrdThresh}}) is infinite, since two points \emph{may}
  interact at any distance apart.

  The function \code{reach(x)} is generic, with methods
  for the case where \code{x} is 
  \itemize{
    \item
    a fitted point process model
    (object of class \code{"ppm"}, usually obtained from the model-fitting
    function \code{\link{ppm}});
    \item
    an interpoint interaction structure (object of class
    \code{"interact"}), created by one of the functions
    \code{\link{Poisson}},
    \code{\link{Strauss}},
    \code{\link{StraussHard}},
    \code{\link{MultiStrauss}},
    \code{\link{MultiStraussHard}},
    \code{\link{Softcore}},
    \code{\link{DiggleGratton}},
    \code{\link{Pairwise}},
    \code{\link{PairPiece}},
    \code{\link{Geyer}},
    \code{\link{LennardJones}},
    \code{\link{Saturated}},
    \code{\link{OrdThresh}}
    or
    \code{\link{Ord}};
    \item
    a point process model for simulation (object of class
    \code{"rmhmodel"}), usually obtained from \code{\link{rmhmodel}}.
  }
  When \code{x} is an \code{"interact"} object,
  \code{reach(x)} returns the maximum possible interaction range
  for any point process model with interaction structure given by \code{x}.
  For example, \code{reach(Strauss(0.2))} returns \code{0.2}.
  
  When \code{x} is a \code{"ppm"} object,
  \code{reach(x)} returns the interaction range
  for the point process model represented by \code{x}.
  For example, a fitted Strauss process model
  with parameters \code{beta,gamma,r} will return
  either \code{0} or \code{r}, depending on whether the fitted
  interaction parameter \code{gamma} is equal or not equal to 1.

  For some point process models, such as the soft core process
  (see \code{\link{Softcore}}), the interaction distance is
  infinite, because the interaction terms are positive for all
  pairs of points. A practical solution is to compute 
  the distance at which the interaction contribution
  from a pair of points falls below a threshold \code{epsilon},
  on the scale of the log conditional intensity. This is done
  by setting the argument \code{epsilon} to a positive value.
}
\seealso{
  \code{\link{ppm}},
    \code{\link{Poisson}},
    \code{\link{Strauss}},
    \code{\link{StraussHard}},
    \code{\link{MultiStrauss}},
    \code{\link{MultiStraussHard}},
    \code{\link{Softcore}},
    \code{\link{DiggleGratton}},
    \code{\link{Pairwise}},
    \code{\link{PairPiece}},
    \code{\link{Geyer}},
    \code{\link{LennardJones}},
    \code{\link{Saturated}},
    \code{\link{OrdThresh}},
    \code{\link{Ord}},
    \code{\link{rmhmodel}}
}
\examples{
    reach(Poisson())
    # returns 0

    reach(Strauss(r=7))
    # returns 7
    data(swedishpines)
    fit <- ppm(swedishpines, ~1, Strauss(r=7))
    reach(fit)
    # returns 7

    reach(OrdThresh(42))
    # returns Inf
    
    reach(MultiStrauss(1:2, matrix(c(1,3,3,1),2,2)))
    # returns 3
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{models}
