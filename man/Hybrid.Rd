\name{Hybrid}
\alias{Hybrid}
\title{
  Hybrid Interaction Point Process Model
}
\description{
  Creates an instance of a hybrid point process model
  which can then be fitted to point pattern data.
}
\usage{
Hybrid(...)
}
\arguments{
  \item{\dots}{
    Two or more interactions (objects of class \code{"interact"})
    or objects which can be converted to interactions.
    See Details.
  }
}
\details{
  A \emph{hybrid} (Baddeley, Turner, Mateu and Bevan, 2013)
  is a point process model created by combining two or more
  point process models, or an interpoint interaction created by combining
  two or more interpoint interactions.
  
  The \emph{hybrid} of two point processes, with probability densities
  \eqn{f(x)} and \eqn{g(x)} respectively,
  is the point process with probability density 
  \deqn{h(x) = c \, f(x) \, g(x)}{h(x) = c * f(x) * g(x)}
  where \eqn{c} is a normalising constant.

  Equivalently, the hybrid of two point processes with conditional intensities
  \eqn{\lambda(u,x)}{lambda(u,x)} and \eqn{\kappa(u,x)}{kappa(u,x)}
  is the point process with conditional intensity
  \deqn{
    \phi(u,x) = \lambda(u,x) \, \kappa(u,x).
  }{
    phi(u,x) = lambda(u,x) * kappa(u,x).
  }
  The hybrid of \eqn{m > 3} point processes is defined in a similar way.
  
  The function \code{\link{ppm}}, which fits point process models to 
  point pattern data, requires an argument 
  of class \code{"interact"} describing the interpoint interaction
  structure of the model to be fitted. 
  The appropriate description of a hybrid interaction is
  yielded by the function \code{Hybrid()}. 
  
  The arguments \code{\dots} will be interpreted as interpoint interactions
  (objects of class \code{"interact"}) and the result will be the hybrid
  of these interactions. Each argument must either be an
  interpoint interaction (object of class \code{"interact"}),
  or a point process model (object of class \code{"ppm"}) from which the
  interpoint interaction will be extracted.

  The arguments \code{\dots} may also be given in the form
  \code{name=value}. This is purely cosmetic: it can be used to attach
  simple mnemonic names to the component interactions, and makes the
  printed output from \code{\link{print.ppm}} neater.
}
\value{
  An object of class \code{"interact"}
  describing an interpoint interaction structure.
}
\references{
  Baddeley, A., Turner, R., Mateu, J. and Bevan, A. (2013)
  Hybrids of Gibbs point process models and their implementation.
  To appear in \emph{Journal of Statistical Software}.
}
\seealso{
  \code{\link{ppm}}
}
\examples{
  Hybrid(Strauss(0.1), Geyer(0.2, 3))

  Hybrid(Ha=Hardcore(0.05), St=Strauss(0.1), Ge=Geyer(0.2, 3))

  fit <- ppm(redwood, ~1, Hybrid(A=Strauss(0.02), B=Geyer(0.1, 2)))
  fit

  ctr <- rmhcontrol(nrep=5e4, expand=1)
  plot(simulate(fit, control=ctr))

  # hybrid components can be models (including hybrid models)
  Hybrid(fit, S=Softcore(0.5))

  # plot.fii only works if every component is a pairwise interaction
  data(swedishpines)
  fit2 <- ppm(swedishpines, ~1, Hybrid(DG=DiggleGratton(2,10), S=Strauss(5)))
  plot(fitin(fit2))
  plot(fitin(fit2), separate=TRUE, mar.panel=rep(4,4))
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{models}
