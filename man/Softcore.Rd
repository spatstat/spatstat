\name{Softcore}
\alias{Softcore}
\title{The Soft Core Point Process Model}
\description{
Creates an instance of the Soft Core point process model
which can then be fitted to point pattern data.
}
\usage{
  Softcore(kappa, sigma0=NA)
}
\arguments{
  \item{kappa}{The exponent \eqn{\kappa}{kappa} of the Soft Core
    interaction}
  \item{sigma0}{
    Optional. Initial estimate of the parameter \eqn{\sigma}{sigma}.
    A positive number.
  }
}
\value{
  An object of class \code{"interact"}
  describing the interpoint interaction
  structure of the Soft Core process with exponent \eqn{\kappa}{kappa}.
}
\details{
  The (stationary)
  Soft Core point process with parameters \eqn{\beta}{beta} and
  \eqn{\sigma}{sigma} and exponent \eqn{\kappa}{kappa}
  is the pairwise interaction point process in which 
  each point contributes a factor \eqn{\beta}{beta} to the 
  probability density of the point pattern, and each pair of points
  contributes a factor
  \deqn{
    \exp \left\{ - \left( \frac{\sigma}{d} \right)^{2/\kappa} \right\}
  }{
    exp( - (sigma/d)^(2/kappa) )
  }
  to the density, where \eqn{d} is the distance between the two points.

  Thus the process has probability density
  \deqn{
    f(x_1,\ldots,x_n) =
    \alpha \beta^{n(x)}
    \exp \left\{ - \sum_{i < j} \left(
                 \frac{\sigma}{||x_i-x_j||}
    \right)^{2/\kappa} \right\}
  }{
    f(x_1,\ldots,x_n) =
    alpha . beta^n(x) exp( - sum (sigma/||x[i]-x[j]||)^(2/kappa))
  }
  where \eqn{x_1,\ldots,x_n}{x[1],\ldots,x[n]} represent the 
  points of the pattern, \eqn{n(x)} is the number of points in the
  pattern, \eqn{\alpha}{alpha} is the normalising constant,
  and the sum on the right hand side is
  over all unordered pairs of points of the pattern. 

  This model describes an ``ordered'' or ``inhibitive'' process,
  with the interpoint interaction decreasing smoothly with distance.
  The strength of interaction is controlled by the
  parameter \eqn{\sigma}{sigma}, a positive real number,
  with larger values corresponding
  to stronger interaction; and by the exponent \eqn{\kappa}{kappa}
  in the range \eqn{(0,1)}, with larger values corresponding to
  weaker interaction.
  If \eqn{\sigma = 0}{sigma = 0}
  the model reduces to the Poisson point process.
  If \eqn{\sigma > 0}{sigma > 0},
  the process is well-defined only for \eqn{\kappa}{kappa} in \eqn{(0,1)}.
  The limit of the model as \eqn{\kappa \to 0}{kappa -> 0} is the
  hard core process with hard core distance \eqn{h=\sigma}{h=sigma}.
 
  The nonstationary Soft Core process is similar except that 
  the contribution of each individual point \eqn{x_i}{x[i]}
  is a function \eqn{\beta(x_i)}{beta(x[i])}
  of location, rather than a constant beta. 
 
  The function \code{\link{ppm}()}, which fits point process models to 
  point pattern data, requires an argument 
  of class \code{"interact"} describing the interpoint interaction
  structure of the model to be fitted. 
  The appropriate description of the Soft Core process pairwise interaction is
  yielded by the function \code{Softcore()}. See the examples below.
 
  The main argument is the exponent  \code{kappa}.
  When \code{kappa} is fixed, the model becomes an exponential family
  with canonical parameters \eqn{\log \beta}{log(beta)}
  and \deqn{
    \log \gamma = \frac{2}{\kappa} \log\sigma
  }{
    log(gamma) = (2/kappa) log(sigma)
  }
  The canonical parameters are estimated by \code{\link{ppm}()}, not fixed in
  \code{Softcore()}. 

  The optional argument \code{sigma0} can be used to improve
  numerical stability. If \code{sigma0} is given, it should be a positive
  number, and it should be a rough estimate of the
  parameter \eqn{\sigma}{sigma}. 
}
\seealso{
  \code{\link{ppm}},
  \code{\link{pairwise.family}},
  \code{\link{ppm.object}}
}
\references{
  Ogata, Y, and Tanemura, M. (1981).
  Estimation of interaction potentials of spatial point patterns
  through the maximum likelihood procedure.
  \emph{Annals of the Institute of Statistical Mathematics}, B
  \bold{33}, 315--338.

  Ogata, Y, and Tanemura, M. (1984).
  Likelihood analysis of spatial point patterns.
  \emph{Journal of the Royal Statistical Society, series B}
  \bold{46}, 496--518.
}
\examples{
   data(cells) 
   ppm(cells, ~1, Softcore(kappa=0.5), correction="isotropic")
   # fit the stationary Soft Core process to `cells'
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{models}
