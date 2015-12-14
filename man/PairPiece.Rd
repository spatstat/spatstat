\name{PairPiece}
\alias{PairPiece}
\title{The Piecewise Constant Pairwise Interaction Point Process Model}
\description{
  Creates an instance of a pairwise interaction point process model
  with piecewise constant potential function. The model
  can then be fitted to point pattern data.
}
\usage{
  PairPiece(r)
}
\arguments{
  \item{r}{vector of jump points for the potential function}
}
\value{
  An object of class \code{"interact"}
  describing the interpoint interaction
  structure of a point process. The process is a pairwise interaction process,
  whose interaction potential is piecewise constant, with jumps
  at the distances given in the vector \eqn{r}.
}
\details{
  A pairwise interaction point process in a bounded region
  is a stochastic point process with probability density of the form
  \deqn{
    f(x_1,\ldots,x_n) =
    \alpha \prod_i b(x_i) \prod_{i < j} h(x_i, x_j)
  }{
    f(x_1,\ldots,x_n) =
    alpha . product { b(x[i]) } product { h(x_i, x_j) }
  }
  where \eqn{x_1,\ldots,x_n}{x[1],\ldots,x[n]} represent the 
  points of the pattern. The first product on the right hand side is
  over all points of the pattern; the second product is over all
  unordered pairs of points of the pattern.

  Thus each point \eqn{x_i}{x[i]} of the pattern contributes a factor 
  \eqn{b(x_i)}{b(x[i])} to the probability density, and each pair of
  points \eqn{x_i, x_j}{x[i], x[j]} contributes a factor
  \eqn{h(x_i,x_j)}{h(x[i], x[j])} to the density.

  The pairwise interaction term \eqn{h(u, v)} is called
  \emph{piecewise constant}
  if it depends only on the distance between \eqn{u} and \eqn{v},
  say \eqn{h(u,v) = H(||u-v||)}, and \eqn{H} is a piecewise constant
  function (a function which is constant except for jumps at a finite
  number of places). The use of piecewise constant interaction terms
  was first suggested by Takacs (1986).
 
  The function \code{\link{ppm}()}, which fits point process models to 
  point pattern data, requires an argument 
  of class \code{"interact"} describing the interpoint interaction
  structure of the model to be fitted. 
  The appropriate description of the piecewise constant pairwise
  interaction is yielded by the function \code{PairPiece()}.
  See the examples below.

  The entries of \code{r} must be strictly increasing, positive numbers.
  They are interpreted as the points of discontinuity of \eqn{H}.
  It is assumed that \eqn{H(s) =1} for all \eqn{s > r_{max}}{s > rmax}
  where \eqn{r_{max}}{rmax} is the maximum value in \code{r}. Thus the
  model has as many regular parameters (see \code{\link{ppm}}) 
  as there are entries in \code{r}. The \eqn{i}-th regular parameter
  \eqn{\theta_i}{theta[i]} is the logarithm of the value of the
  interaction function \eqn{H} on the interval
  \eqn{[r_{i-1},r_i)}{[r[i-1],r[i])}. 

  If \code{r} is a single number, this model is similar to the 
  Strauss process, see \code{\link{Strauss}}. The difference is that
  in \code{PairPiece} the interaction function is continuous on the
  right, while in \code{\link{Strauss}} it is continuous on the left.

  The analogue of this model for multitype point processes
  has not yet been implemented.
}
\seealso{
  \code{\link{ppm}},
  \code{\link{pairwise.family}},
  \code{\link{ppm.object}},
  \code{\link{Strauss}}
  \code{\link{rmh.ppm}}
}
\examples{
   PairPiece(c(0.1,0.2))
   # prints a sensible description of itself
   data(cells) 

   \dontrun{
   ppm(cells, ~1, PairPiece(r = c(0.05, 0.1, 0.2)))
   # fit a stationary piecewise constant pairwise interaction process
   }

   ppm(cells, ~polynom(x,y,3), PairPiece(c(0.05, 0.1)))
   # nonstationary process with log-cubic polynomial trend
}
\references{
  Takacs, R. (1986)
  Estimator for the pair potential of a Gibbsian point process.
  \emph{Statistics} \bold{17}, 429--433.
}
\author{Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  and Rolf Turner \email{r.turner@auckland.ac.nz}
  
}
\keyword{spatial}
\keyword{models}
