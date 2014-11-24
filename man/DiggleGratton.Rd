\name{DiggleGratton}
\alias{DiggleGratton}
\title{Diggle-Gratton model}
\description{
  Creates an instance of the Diggle-Gratton pairwise interaction
  point process model, which can then be fitted to point pattern data.
}
\usage{
  DiggleGratton(delta=NA, rho)
}
\arguments{
  \item{delta}{lower threshold \eqn{\delta}{delta}}
  \item{rho}{upper threshold \eqn{\rho}{rho}}
}
\value{
  An object of class \code{"interact"}
  describing the interpoint interaction
  structure of a point process. 
}
\details{
  Diggle and Gratton (1984, pages 208-210)
  introduced the pairwise interaction point
  process with pair potential \eqn{h(t)} of the form
  \deqn{
    h(t) = \left( \frac{t-\delta}{\rho-\delta} \right)^\kappa
    \quad\quad \mbox{  if  } \delta \le t \le \rho
  }{
    h(t) = ((t - delta)/(rho - delta))^kappa, {    } delta <= t <= rho
  }
  with \eqn{h(t) = 0} for \eqn{t < \delta}{t < delta}
  and  \eqn{h(t) = 1} for \eqn{t > \rho}{t > rho}.
  Here \eqn{\delta}{delta}, \eqn{\rho}{rho} and \eqn{\kappa}{kappa}
  are parameters.

  Note that we use the symbol \eqn{\kappa}{kappa}
  where Diggle and Gratton (1984) and Diggle, Gates and Stibbard (1987)
  use \eqn{\beta}{beta}, since in \pkg{spatstat} we reserve the symbol
  \eqn{\beta}{beta} for an intensity parameter.

  The parameters must all be nonnegative,
  and must satisfy \eqn{\delta \le \rho}{delta <= rho}.

  The potential is inhibitory, i.e.\ this model is only appropriate for
  regular point patterns. The strength of inhibition increases with
  \eqn{\kappa}{kappa}. For \eqn{\kappa=0}{kappa=0} the model is
  a hard core process with hard core radius \eqn{\delta}{delta}.
  For \eqn{\kappa=\infty}{kappa=infinity} the model is a hard core
  process with hard core radius \eqn{\rho}{rho}.

  The irregular parameters 
  \eqn{\delta, \rho}{delta, rho} must be given in the call to
  \code{DiggleGratton}, while the
  regular parameter \eqn{\kappa}{kappa} will be estimated.

  If the lower threshold \code{delta} is missing or \code{NA},
  it will be estimated from the data when \code{\link{ppm}} is called.
  The estimated value of \code{delta} is the minimum nearest neighbour distance
  multiplied by \eqn{n/(n+1)}, where \eqn{n} is the
  number of data points.
}
\seealso{
  \code{\link{ppm}},
  \code{\link{ppm.object}},
  \code{\link{Pairwise}}
}
\examples{
   ppm(cells ~1, DiggleGratton(0.05, 0.1))
}
\references{
  Diggle, P.J., Gates, D.J. and Stibbard, A. (1987)
  A nonparametric estimator for pairwise-interaction point processes.
  \emph{Biometrika} \bold{74}, 763 -- 770.

  Diggle, P.J. and Gratton, R.J. (1984)
  Monte Carlo methods of inference for implicit statistical models.
  \emph{Journal of the Royal Statistical Society, series B}
  \bold{46}, 193 -- 212.
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{models}

