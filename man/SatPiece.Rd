\name{SatPiece}
\alias{SatPiece}
\title{Piecewise Constant Saturated Pairwise Interaction Point Process Model}
\description{
  Creates an instance of a saturated pairwise interaction point process model
  with piecewise constant potential function. The model
  can then be fitted to point pattern data.
}
\usage{
  SatPiece(r, sat)
}
\arguments{
  \item{r}{vector of jump points for the potential function}
  \item{sat}{
    vector of saturation values,
    or a single saturation value
  }
}
\value{
  An object of class \code{"interact"}
  describing the interpoint interaction
  structure of a point process. 
}
\details{
  This is a generalisation of the Geyer saturation point process model,
  described in \code{\link{Geyer}}, to the case of multiple interaction
  distances. It can also be described as the saturated analogue of a
  pairwise interaction process with piecewise-constant pair potential,
  described in \code{\link{PairPiece}}.

  The saturated point process with interaction radii
  \eqn{r_1,\ldots,r_k}{r[1], ..., r[k]},
  saturation thresholds \eqn{s_1,\ldots,s_k}{s[1],...,s[k]},
  intensity parameter \eqn{\beta}{beta} and
  interaction parameters
  \eqn{\gamma_1,\ldots,gamma_k}{gamma[1], ..., gamma[k]},
  is the point process
  in which each point
  \eqn{x_i}{x[i]} in the pattern \eqn{X}
  contributes a factor
  \deqn{
    \beta \gamma_1^{v_1(x_i, X)} \ldots gamma_k^{v_k(x_i,X)}
  }{
    beta gamma[1]^v(1, x_i, X) ... gamma[k]^v(k, x_i, X)
  }
  to the probability density of the point pattern,
  where
  \deqn{
    v_j(x_i, X) = \min( s_j, t_j(x_i,X) )
  }{
    v(j, x_i, X) = min(s[j], t(j, x_i, X))
  }
  where \eqn{t_j(x_i, X)}{t(j,x[i],X)} denotes the
  number of points in the pattern \eqn{X} which lie
  at a distance between \eqn{r_{j-1}}{r[j-1]} and \eqn{r_j}{r[j]}
  from the point \eqn{x_i}{x[i]}. We take \eqn{r_0 = 0}{r[0] = 0}
  so that \eqn{t_1(x_i,X)}{t(1, x[i], X)} is the number of points of
  \eqn{X} that lie within a distance \eqn{r_1}{r[1]} of the point
  \eqn{x_i}{x[i]}. 

  \code{SatPiece} is used to fit this model to data.
  The function \code{\link{ppm}()}, which fits point process models to 
  point pattern data, requires an argument 
  of class \code{"interact"} describing the interpoint interaction
  structure of the model to be fitted. 
  The appropriate description of the piecewise constant Saturated pairwise
  interaction is yielded by the function \code{SatPiece()}.
  See the examples below.

  Simulation of this point process model is not yet implemented.
  This model is not locally stable (the conditional intensity is
  unbounded).
  
  The argument \code{r} specifies the vector of interaction distances.
  The entries of \code{r} must be strictly increasing, positive numbers.

  The argument \code{sat} specifies the vector of saturation parameters.
  It should be a vector of the same length as \code{r}, and its entries
  should be nonnegative numbers. Thus \code{sat[1]} corresponds to the
  distance range from \code{0} to \code{r[1]}, and \code{sat[2]} to the
  distance range from \code{r[1]} to \code{r[2]}, etc.
  Alternatively \code{sat} may be a single number, and this saturation
  value will be applied to every distance range.

  Infinite values of the
  saturation parameters are also permitted; in this case
  \eqn{v_j(x_i,X) = t_j(x_i,X)}{v(j, x_i, X) = t(j, x_i, X)}
  and there is effectively no `saturation' for the distance range in
  question. If all the saturation parameters are set to \code{Inf} then
  the model is effectively a pairwise interaction process, equivalent to
  \code{\link{PairPiece}} (however the interaction parameters
  \eqn{\gamma}{gamma} obtained from \code{\link{SatPiece}} are the
  square roots of the parameters \eqn{\gamma}{gamma}
  obtained from \code{\link{PairPiece}}).
   
  If \code{r} is a single number, this model is virtually equivalent to the 
  Geyer process, see \code{\link{Geyer}}. 
}
\seealso{
  \code{\link{ppm}},
  \code{\link{pairsat.family}},
  \code{\link{Geyer}},
  \code{\link{PairPiece}},
  \code{\link{BadGey}}.
}
\examples{
   SatPiece(c(0.1,0.2), c(1,1))
   # prints a sensible description of itself
   SatPiece(c(0.1,0.2), 1)
   data(cells) 
   ppm(cells, ~1, SatPiece(c(0.07, 0.1, 0.13), 2))
   # fit a stationary piecewise constant Saturated pairwise interaction process

   \dontrun{
   ppm(cells, ~polynom(x,y,3), SatPiece(c(0.07, 0.1, 0.13), 2))
   # nonstationary process with log-cubic polynomial trend
   }
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
  in collaboration with Hao Wang and Jeff Picka
}
\keyword{spatial}
\keyword{models}
