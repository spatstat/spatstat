\name{Fiksel}
\alias{Fiksel}
\title{The Fiksel Interaction}
\description{
  Creates an instance of Fiksel's double exponential
  pairwise interaction point process model,
  which can then be fitted to point pattern data.
}
\usage{
  Fiksel(r, hc=NA, kappa)
}
\arguments{
  \item{r}{The interaction radius of the Fiksel model}
  \item{hc}{The hard core distance}
  \item{kappa}{The rate parameter}
}
\value{
  An object of class \code{"interact"}
  describing the interpoint interaction
  structure of the Fiksel
  process with interaction radius \eqn{r},
  hard core distance \code{hc} and
  rate parameter \code{kappa}.
}
\details{
  Fiksel (1984) introduced a pairwise interaction point process
  with the following interaction function \eqn{c}.
  For two points \eqn{u} and \eqn{v} separated by a distance
  \eqn{d=||u-v||}, the interaction
  \eqn{c(u,v)} is equal to \eqn{0} if \eqn{d < h},
  equal to \eqn{1} if \eqn{d > r}, and
  equal to
  \deqn{ \exp(a \exp(-\kappa d))}{exp(a * exp(-kappa * d))}
  if \eqn{h \le d \le r}{h <= d <= r}, where
  \eqn{h,r,\kappa,a}{h,r,kappa,a} are parameters.
  
  A graph of this interaction function is shown in the Examples.
  The interpretation of the parameters is as follows.
  \itemize{
    \item \eqn{h} is the hard core distance: distinct points are
    not permitted to come closer than a distance \eqn{h} apart.
    \item \eqn{r} is the interaction range: points further than
    this distance do not interact.
    \item \eqn{\kappa}{kappa} is the rate or slope parameter,
    controlling the decay of the interaction as distance increases.
    \item \eqn{a} is the interaction strength parameter,
    controlling the strength and type of interaction.
    If \eqn{a} is zero, the process is Poisson. If \code{a} is positive,
    the process is clustered. If \code{a} is negative, the process is
    inhibited (regular).
  }
  
  The function \code{\link{ppm}()}, which fits point process models to 
  point pattern data, requires an argument 
  of class \code{"interact"} describing the interpoint interaction
  structure of the model to be fitted. 
  The appropriate description of the Fiksel
  pairwise interaction is
  yielded by the function \code{Fiksel()}. See the examples below.
  
  The parameters \eqn{h}, \eqn{r} and \eqn{\kappa}{kappa} must be
  fixed and given in the call to \code{Fiksel}, while the canonical
  parameter \eqn{a} is estimated by \code{\link{ppm}()}.
  
  To estimate \eqn{h}, \eqn{r} and\eqn{\kappa}{kappa}
  it is possible to use \code{\link{profilepl}}. The maximum likelihood
  estimator of\eqn{h} is the minimum interpoint distance.

  If the hard core distance argument \code{hc} is missing or \code{NA},
  it will be estimated from the data when \code{\link{ppm}} is called.
  The estimated value of \code{hc} is the minimum nearest neighbour distance
  multiplied by \eqn{n/(n+1)}, where \eqn{n} is the
  number of data points.
  
  See also Stoyan, Kendall and Mecke (1987) page 161.
}
\seealso{
  \code{\link{ppm}},
  \code{\link{pairwise.family}},
  \code{\link{ppm.object}},
  \code{\link{StraussHard}}
}
\references{
  Baddeley, A. and Turner, R. (2000)
  Practical maximum pseudolikelihood for spatial point patterns.
  \emph{Australian and New Zealand Journal of Statistics}
  \bold{42}, 283--322.

Fiksel, T. (1984)
Estimation of parameterized pair potentials
of marked and non-marked Gibbsian point processes.
\emph{Electronische Informationsverabeitung und Kybernetika}
\bold{20}, 270--278.

 Stoyan, D, Kendall, W.S. and Mecke, J. (1987)
 \emph{Stochastic geometry and its applications}.  Wiley.
}
\examples{
   Fiksel(r=1,hc=0.02, kappa=2)
   # prints a sensible description of itself

   data(spruces)
   X <- unmark(spruces)

   fit <- ppm(X ~ 1, Fiksel(r=3.5, kappa=1))
   plot(fitin(fit))
}
\author{
  Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  and
  Rolf Turner \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{models}
