\name{Strauss}
\alias{Strauss}
\title{The Strauss Point Process Model}
\description{
Creates an instance of the Strauss point process model
which can then be fitted to point pattern data.
}
\usage{
  Strauss(r)
}
\arguments{
  \item{r}{The interaction radius of the Strauss process}
}
\value{
  An object of class \code{"interact"}
  describing the interpoint interaction
  structure of the Strauss process with interaction radius \eqn{r}.
}
\details{
  The (stationary) Strauss process with interaction radius \eqn{r} and 
  parameters \eqn{\beta}{beta} and \eqn{\gamma}{gamma}
  is the pairwise interaction point process
  in which each point contributes a factor \eqn{\beta}{beta} to the 
  probability density of the point pattern, and each pair of points
  closer than \eqn{r} units apart contributes a factor
  \eqn{\gamma}{gamma} to the density.

  Thus the probability density is
  \deqn{
    f(x_1,\ldots,x_n) =
    \alpha \beta^{n(x)} \gamma^{s(x)}
  }{
    f(x_1,\ldots,x_n) =
    alpha . beta^n(x) gamma^s(x)
  }
  where \eqn{x_1,\ldots,x_n}{x[1],\ldots,x[n]} represent the 
  points of the pattern, \eqn{n(x)} is the number of points in the
  pattern, \eqn{s(x)} is the number of distinct unordered pairs of
  points that are closer than \eqn{r} units apart,
  and \eqn{\alpha}{alpha} is the normalising constant.

  The interaction parameter \eqn{\gamma}{gamma} must be less than
  or equal to \eqn{1}
  so that this model describes an ``ordered'' or ``inhibitive'' pattern.
 
  The nonstationary Strauss process is similar except that 
  the contribution of each individual point \eqn{x_i}{x[i]}
  is a function \eqn{\beta(x_i)}{beta(x[i])}
  of location, rather than a constant beta. 
 
  The function \code{\link{ppm}()}, which fits point process models to 
  point pattern data, requires an argument 
  of class \code{"interact"} describing the interpoint interaction
  structure of the model to be fitted. 
  The appropriate description of the Strauss process pairwise interaction is
  yielded by the function \code{Strauss()}. See the examples below.
 
  Note the only argument is the interaction radius \code{r}.
  When \code{r} is fixed, the model becomes an exponential family.
  The canonical parameters \eqn{\log(\beta)}{log(beta)}
  and \eqn{\log(\gamma)}{log(gamma)}
  are estimated by \code{\link{ppm}()}, not fixed in
  \code{Strauss()}. 
}
\seealso{
  \code{\link{ppm}},
  \code{\link{pairwise.family}},
  \code{\link{ppm.object}}
}
\references{
 Kelly, F.P. and Ripley, B.D. (1976)
 On Strauss's model for clustering.
 \emph{Biometrika} \bold{63}, 357--360.

 Strauss, D.J. (1975)
 A model for clustering.
 \emph{Biometrika} \bold{63}, 467--475.
}
\examples{
   Strauss(r=0.1)
   # prints a sensible description of itself
   data(cells)

   \dontrun{
   ppm(cells, ~1, Strauss(r=0.07))
   # fit the stationary Strauss process to `cells'
   }


   ppm(cells, ~polynom(x,y,3), Strauss(r=0.07))
   # fit a nonstationary Strauss process with log-cubic polynomial trend
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{models}
