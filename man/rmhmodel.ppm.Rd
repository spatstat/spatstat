\name{rmhmodel.ppm}
\alias{rmhmodel.ppm}
\title{Interpret Fitted Model for Metropolis-Hastings Simulation.}
\description{
  Converts a fitted point process model
  into a format that can be used to simulate the model
  by the Metropolis-Hastings algorithm. 
}
\usage{
  \method{rmhmodel}{ppm}(model, win, ..., verbose=TRUE, project=TRUE,
                         control=rmhcontrol(),
                         new.coef=NULL)
}
\arguments{
  \item{model}{
    Fitted point process model (object of class \code{"ppm"}).
  }
  \item{win}{
    Optional. Window in which the simulations should be generated.
  }
  \item{\dots}{Ignored.}
  \item{verbose}{
    Logical flag indicating whether to print progress reports
    while the model is being converted.
  }
  \item{project}{Logical flag indicating what to do if the fitted model
    does not correspond to a valid point process. See Details.}
  \item{control}{
    Parameters determining the iterative behaviour of the simulation
    algorithm. Passed to \code{\link{rmhcontrol}}.
  }
  \item{new.coef}{
    New values for the canonical parameters of the model.
    A numeric vector of the same length as \code{coef(model)}.
  }
}
\value{
  An object of class \code{"rmhmodel"}, which is essentially
  a list of parameter values for the model.
  
  There is a \code{print} method for this class, which prints
  a sensible description of the model chosen.
}
\details{
  The generic function \code{\link{rmhmodel}} takes a
  description of a point process model in some format, and
  converts it into an object of class \code{"rmhmodel"}
  so that simulations of the model can be generated using
  the Metropolis-Hastings algorithm \code{\link{rmh}}. 
  
  This function \code{rmhmodel.ppm} is the method for
  the class \code{"ppm"} of fitted point process models.

  The argument \code{model} should be a fitted point process model
  (object of class \code{"ppm"}) typically obtained from the
  model-fitting function \code{\link{ppm}}.
  This will be converted into an object of class \code{"rmhmodel"}.

  The optional argument \code{win} specifies the window
  in which the pattern is to be generated.  If specified, it must be in
  a form which can be coerced to an object of class \code{owin}
  by \code{\link{as.owin}}.

  Not all fitted point process models
  obtained from \code{\link{ppm}} can be simulated.
  We have not yet implemented simulation code for
  the \code{\link{LennardJones}} and \code{\link{OrdThresh}}
  models.

  It is also possible that a fitted point process model
  obtained from \code{\link{ppm}}  may not correspond to a valid
  point process. For example a fitted model with the \code{\link{Strauss}}
  interpoint interaction may have any value of the interaction parameter
  \eqn{\gamma}{gamma}; however the Strauss 
  process is not well-defined for 
  \eqn{\gamma > 1}{gamma > 1} (Kelly and Ripley, 1976).

  The argument \code{project} determines what to do in such cases.
  If \code{project=FALSE}, a fatal error will occur.
  If \code{project=TRUE}, the fitted model parameters will be
  adjusted to the nearest values which do correspond to a valid
  point process. For example a Strauss process with \eqn{\gamma >
    1}{gamma > 1} will be projected to a Strauss process with
  \eqn{\gamma = 1}{gamma = 1}, equivalent to a Poisson process.
}  
\references{
   Diggle, P. J. (2003) \emph{Statistical Analysis of Spatial Point
   Patterns} (2nd ed.) Arnold, London.

   Diggle, P.J. and Gratton, R.J. (1984)
   Monte Carlo methods of inference for implicit statistical models.
   \emph{Journal of the Royal Statistical Society, series B}
   \bold{46}, 193 -- 212.

   Geyer, C.J. (1999)
   Likelihood Inference for Spatial Point
   Processes. Chapter 3 in  O.E. Barndorff-Nielsen, W.S. Kendall and
   M.N.M. Van Lieshout (eds) \emph{Stochastic Geometry: Likelihood and
   Computation}, Chapman and Hall / CRC,  Monographs on Statistics and
   Applied Probability, number 80. Pages 79--140.

   Kelly, F.P. and Ripley, B.D. (1976)
   On Strauss's model for clustering.
   \emph{Biometrika} \bold{63}, 357--360.
}
\seealso{
  \code{\link{rmhmodel}},
  \code{\link{rmhmodel.list}},
  \code{\link{rmhmodel.default}},
  \code{\link{rmh}},
  \code{\link{rmhcontrol}},
  \code{\link{rmhstart}},
  \code{\link{ppm}},
  \code{\link{AreaInter}},
  \code{\link{BadGey}},
  \code{\link{DiggleGatesStibbard}},
  \code{\link{DiggleGratton}},
  \code{\link{Fiksel}},
  \code{\link{Geyer}},
  \code{\link{Hardcore}},
  \code{\link{Hybrid}},
  \code{\link{LennardJones}},
  \code{\link{MultiStrauss}},
  \code{\link{MultiStraussHard}},
  \code{\link{PairPiece}},
  \code{\link{Poisson}},
  \code{\link{Softcore}},
  \code{\link{Strauss}},
  \code{\link{StraussHard}},
  \code{\link{Triplets}}
}
\examples{
  data(cells)
  fit1 <- ppm(cells, ~1, Strauss(0.07))
  mod1 <- rmhmodel(fit1)

  fit2 <- ppm(cells, ~x, Geyer(0.07, 2))
  mod2 <- rmhmodel(fit2)

  fit3 <- ppm(cells, ~x, Hardcore(0.07))
  mod3 <- rmhmodel(fit3)

  # Then rmh(mod1), etc
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{datagen}

