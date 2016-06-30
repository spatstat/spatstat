\name{DiggleGatesStibbard}
\alias{DiggleGatesStibbard}
\title{Diggle-Gates-Stibbard Point Process Model}
\description{
Creates an instance of the Diggle-Gates-Stibbard point process model
which can then be fitted to point pattern data.
}
\usage{
  DiggleGatesStibbard(rho)
}
\arguments{
  \item{rho}{Interaction range}
}
\value{
  An object of class \code{"interact"}
  describing the interpoint interaction
  structure of the Diggle-Gates-Stibbard
  process with interaction range \code{rho}.
}
\details{
  Diggle, Gates and Stibbard (1987) proposed a pairwise interaction
  point process in which each pair of points separated by
  a distance \eqn{d} contributes a factor \eqn{e(d)} to the
  probability density, where
  \deqn{
    e(d) = \sin^2\left(\frac{\pi d}{2\rho}\right)
  }{
    e(d) = sin^2((pi * d)/(2 * rho))
  }
  for \eqn{d < \rho}{d < rho}, and \eqn{e(d)} is equal to 1
  for \eqn{d \ge \rho}{d >= rho}.

  The function \code{\link{ppm}()}, which fits point process models to 
  point pattern data, requires an argument 
  of class \code{"interact"} describing the interpoint interaction
  structure of the model to be fitted. 
  The appropriate description of the Diggle-Gates-Stibbard 
  pairwise interaction is 
  yielded by the function \code{DiggleGatesStibbard()}.
  See the examples below.

  Note that this model does not have any regular parameters
  (as explained in the section on Interaction Parameters
  in the help file for \code{\link{ppm}}).
  The parameter \eqn{rho} is not estimated by \code{\link{ppm}}.
}
\seealso{
  \code{\link{ppm}},
  \code{\link{pairwise.family}},
  \code{\link{DiggleGratton}},
  \code{\link{rDGS}},
  \code{\link{ppm.object}}
}
\references{
  Baddeley, A. and Turner, R. (2000)
  Practical maximum pseudolikelihood for spatial point patterns.
  \emph{Australian and New Zealand Journal of Statistics}
  \bold{42}, 283--322.

  Ripley, B.D. (1981) 
  \emph{Spatial statistics}.
  John Wiley and Sons.
  
  Diggle, P.J., Gates, D.J., and Stibbard, A. (1987)
   A nonparametric estimator for pairwise-interaction point processes.
   Biometrika \bold{74}, 763 -- 770.
   \emph{Scandinavian Journal of Statistics} \bold{21}, 359--373.
}
\examples{
   DiggleGatesStibbard(0.02)
   # prints a sensible description of itself

   \dontrun{
   ppm(cells ~1, DiggleGatesStibbard(0.05))
   # fit the stationary D-G-S process to `cells'
   }

   ppm(cells ~ polynom(x,y,3), DiggleGatesStibbard(0.05))
   # fit a nonstationary D-G-S process
   # with log-cubic polynomial trend
}
\author{
  \spatstatAuthors
}
\keyword{spatial}
\keyword{models}
