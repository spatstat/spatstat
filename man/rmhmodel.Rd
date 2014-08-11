\name{rmhmodel}
\alias{rmhmodel}
\title{Define Point Process Model for Metropolis-Hastings Simulation.}
\description{
  Builds a description of a point process model
  for use in simulating the model by the Metropolis-Hastings
  algorithm. 
}
\usage{
  rmhmodel(...)
}
\arguments{
  \item{\dots}{Arguments specifying the point process model
    in some format.
  }
}
\value{
  An object of class \code{"rmhmodel"}, which is essentially
  a list of parameter values for the model.
  
  There is a \code{print} method for this class, which prints
  a sensible description of the model chosen.
}
\details{
  Simulated realisations of many point process models
  can be generated using the Metropolis-Hastings algorithm
  \code{\link{rmh}}. The algorithm requires the model to be specified
  in a particular format: an object of class \code{"rmhmodel"}.

  The function \code{\link{rmhmodel}} takes a
  description of a point process model in some other format, and
  converts it into an object of class \code{"rmhmodel"}.
  It also checks that the parameters of the model are valid.

  The function \code{\link{rmhmodel}} is generic, with methods
  for
  \describe{
    \item{fitted point process models:}{
      an object of class \code{"ppm"}, obtained by a call to the
      model-fitting function \code{\link{ppm}}.
      See \code{\link{rmhmodel.ppm}}.
    }
    \item{lists:}{
      a list of parameter values in a certain format.
      See \code{\link{rmhmodel.list}}.
    }
    \item{default:}{
      parameter values specified as separate arguments to \code{\dots}.
      See \code{\link{rmhmodel.default}}.
    }
  }
}
\references{
   Diggle, P. J. (2003) \emph{Statistical Analysis of Spatial Point
   Patterns} (2nd ed.) Arnold, London.

   Diggle, P.J. and Gratton, R.J. (1984)
   Monte Carlo methods of inference for implicit statistical models.
   \emph{Journal of the Royal Statistical Society, series B}
   \bold{46}, 193 -- 212.

   Diggle, P.J., Gates, D.J., and Stibbard, A. (1987)
   A nonparametric estimator for pairwise-interaction point processes.
   Biometrika \bold{74}, 763 -- 770.
   \emph{Scandinavian Journal of Statistics} \bold{21}, 359--373.

   Geyer, C.J. (1999)
   Likelihood Inference for Spatial Point
   Processes. Chapter 3 in  O.E. Barndorff-Nielsen, W.S. Kendall and
   M.N.M. Van Lieshout (eds) \emph{Stochastic Geometry: Likelihood and
   Computation}, Chapman and Hall / CRC,  Monographs on Statistics and
   Applied Probability, number 80. Pages 79--140.
}

\seealso{
  \code{\link{rmhmodel.ppm}},
  \code{\link{rmhmodel.default}},
  \code{\link{rmhmodel.list}},
  \code{\link{rmh}},
  \code{\link{rmhcontrol}},
  \code{\link{rmhstart}},
  \code{\link{ppm}},
  \code{\link{Strauss}},
  \code{\link{Softcore}},
  \code{\link{StraussHard}},
  \code{\link{Triplets}},
  \code{\link{MultiStrauss}},
  \code{\link{MultiStraussHard}},
  \code{\link{DiggleGratton}},
  \code{\link{PairPiece}}
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{datagen}

