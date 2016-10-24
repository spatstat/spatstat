\name{BadGey}
\alias{BadGey}
\title{Hybrid Geyer Point Process Model}
\description{
  Creates an instance of the Baddeley-Geyer point process model, defined
  as a hybrid of several Geyer interactions. The model
  can then be fitted to point pattern data.
}
\usage{
  BadGey(r, sat)
}
\arguments{
  \item{r}{vector of interaction radii}
  \item{sat}{
    vector of saturation parameters,
    or a single common value of saturation parameter
  }
}
\value{
  An object of class \code{"interact"}
  describing the interpoint interaction
  structure of a point process. 
}
\details{
  This is Baddeley's generalisation of the
  Geyer saturation point process model,
  described in \code{\link{Geyer}}, to a process with multiple interaction
  distances. 

  The BadGey point process with interaction radii
  \eqn{r_1,\ldots,r_k}{r[1], \ldots, r[k]},
  saturation thresholds \eqn{s_1,\ldots,s_k}{s[1],\ldots,s[k]},
  intensity parameter \eqn{\beta}{\beta} and
  interaction parameters
  \eqn{\gamma_1,\ldots,gamma_k}{\gamma[1], \ldots, \gamma[k]},
  is the point process
  in which each point
  \eqn{x_i}{x[i]} in the pattern \eqn{X}
  contributes a factor
  \deqn{
    \beta \gamma_1^{v_1(x_i, X)} \ldots gamma_k^{v_k(x_i,X)}
  }{
    \beta \gamma[1]^v(1, x_i, X) \ldots \gamma[k]^v(k, x_i, X)
  }
  to the probability density of the point pattern,
  where
  \deqn{
    v_j(x_i, X) = \min( s_j, t_j(x_i,X) )
  }{
    v(j, x[i], X) = min(s[j], t(j, x[i], X))
  }
  where \eqn{t_j(x_i, X)}{t(j,x[i],X)} denotes the
  number of points in the pattern \eqn{X} which lie
  within a distance \eqn{r_j}{r[j]}
  from the point \eqn{x_i}{x[i]}. 

  \code{BadGey} is used to fit this model to data.
  The function \code{\link{ppm}()}, which fits point process models to 
  point pattern data, requires an argument 
  of class \code{"interact"} describing the interpoint interaction
  structure of the model to be fitted. 
  The appropriate description of the piecewise constant Saturated pairwise
  interaction is yielded by the function \code{BadGey()}.
  See the examples below.

  The argument \code{r} specifies the vector of interaction distances.
  The entries of \code{r} must be strictly increasing, positive numbers.

  The argument \code{sat} specifies the vector of saturation parameters
  that are applied to the point counts \eqn{t_j(x_i, X)}{t(j,x[i],X)}.
  It should be a vector of the same length as \code{r}, and its entries
  should be nonnegative numbers. Thus \code{sat[1]} is applied to the
  count of points within a distance \code{r[1]}, and \code{sat[2]} to the
  count of points within a distance \code{r[2]}, etc.
  Alternatively \code{sat} may be a single number, and this saturation
  value will be applied to every count.

  Infinite values of the
  saturation parameters are also permitted; in this case
  \eqn{v_j(x_i,X) = t_j(x_i,X)}{v(j, x[i], X) = t(j, x[i], X)}
  and there is effectively no `saturation' for the distance range in
  question. If all the saturation parameters are set to \code{Inf} then
  the model is effectively a pairwise interaction process, equivalent to
  \code{\link{PairPiece}} (however the interaction parameters
  \eqn{\gamma}{\gamma} obtained from \code{\link{BadGey}}
  have a complicated relationship to the interaction
  parameters \eqn{\gamma}{\gamma} obtained from \code{\link{PairPiece}}).
   
  If \code{r} is a single number, this model is virtually equivalent to the 
  Geyer process, see \code{\link{Geyer}}. 
}
\seealso{
  \code{\link{ppm}},
  \code{\link{pairsat.family}},
  \code{\link{Geyer}},
  \code{\link{PairPiece}},
  \code{\link{SatPiece}}
}
\examples{
   BadGey(c(0.1,0.2), c(1,1))
   # prints a sensible description of itself
   BadGey(c(0.1,0.2), 1)
   data(cells) 

   # fit a stationary Baddeley-Geyer model
   ppm(cells, ~1, BadGey(c(0.07, 0.1, 0.13), 2))

   # nonstationary process with log-cubic polynomial trend
   \dontrun{
   ppm(cells, ~polynom(x,y,3), BadGey(c(0.07, 0.1, 0.13), 2))
   }
}
\author{
  \adrian
  and
  \rolf
  in collaboration with Hao Wang and Jeff Picka
}
\keyword{spatial}
\keyword{models}
