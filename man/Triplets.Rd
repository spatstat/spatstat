\name{Triplets}
\alias{Triplets}
\title{The Triplet Point Process Model}
\description{
Creates an instance of Geyer's triplet interaction point process model
which can then be fitted to point pattern data.
}
\usage{
  Triplets(r)
}
\arguments{
  \item{r}{The interaction radius of the Triplets process}
}
\value{
  An object of class \code{"interact"}
  describing the interpoint interaction
  structure of the Triplets process with interaction radius \eqn{r}.
}
\details{
  The (stationary) Geyer triplet process (Geyer, 1999)
  with interaction radius \eqn{r} and 
  parameters \eqn{\beta}{beta} and \eqn{\gamma}{gamma}
  is the point process
  in which each point contributes a factor \eqn{\beta}{beta} to the 
  probability density of the point pattern, and each triplet of close points
  contributes a factor \eqn{\gamma}{gamma} to the density.
  A triplet of close points is a group of 3 points,
  each pair of which is closer than \eqn{r} units
  apart.

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
  pattern, \eqn{s(x)} is the number of unordered triples of
  points that are closer than \eqn{r} units apart,
  and \eqn{\alpha}{alpha} is the normalising constant.

  The interaction parameter \eqn{\gamma}{gamma} must be less than
  or equal to \eqn{1}
  so that this model describes an ``ordered'' or ``inhibitive'' pattern.
 
  The nonstationary Triplets process is similar except that 
  the contribution of each individual point \eqn{x_i}{x[i]}
  is a function \eqn{\beta(x_i)}{beta(x[i])}
  of location, rather than a constant beta. 
 
  The function \code{\link{ppm}()}, which fits point process models to 
  point pattern data, requires an argument 
  of class \code{"interact"} describing the interpoint interaction
  structure of the model to be fitted. 
  The appropriate description of the Triplets process pairwise interaction is
  yielded by the function \code{Triplets()}. See the examples below.
 
  Note the only argument is the interaction radius \code{r}.
  When \code{r} is fixed, the model becomes an exponential family.
  The canonical parameters \eqn{\log(\beta)}{log(beta)}
  and \eqn{\log(\gamma)}{log(gamma)}
  are estimated by \code{\link{ppm}()}, not fixed in
  \code{Triplets()}.
}
\seealso{
  \code{\link{ppm}},
  \code{\link{triplet.family}},
  \code{\link{ppm.object}}
}
\references{
  Geyer, C.J. (1999)
  Likelihood Inference for Spatial Point Processes.
  Chapter 3 in 
  O.E. Barndorff-Nielsen, W.S. Kendall and M.N.M. Van Lieshout (eds)
  \emph{Stochastic Geometry: Likelihood and Computation},
  Chapman and Hall / CRC, 
  Monographs on Statistics and Applied Probability, number 80.
  Pages 79--140.
}
\examples{
   Triplets(r=0.1)
   # prints a sensible description of itself

   \dontrun{
   ppm(cells, ~1, Triplets(r=0.1))
   # fit the stationary Triplets process to `cells'
   }


   ppm(cells, ~polynom(x,y,3), Triplets(r=0.1))
   # fit a nonstationary Triplets process with log-cubic polynomial trend
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{models}
