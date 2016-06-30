\name{msr}
\alias{msr}
\title{
  Signed or Vector-Valued Measure
}
\description{
  Defines an object representing a signed measure or vector-valued
  measure on a spatial domain.
}
\usage{
  msr(qscheme, discrete, density, check=TRUE)
}
\arguments{
  \item{qscheme}{
    A quadrature scheme (object of class \code{"quad"} usually
    extracted from a fitted point process model).
  }
  \item{discrete}{
    Vector or matrix containing the values (masses) of the discrete component
    of the measure, for each of the data points in \code{qscheme}.
  }
  \item{density}{
    Vector or matrix containing values of the density of the
    diffuse component of the measure, for each of the
    quadrature points in \code{qscheme}.
  }
  \item{check}{
    Logical. Whether to check validity of the arguments.
  }
}
\details{
  This function creates an object that represents a
  signed or vector valued \emph{measure} on the two-dimensional plane.
  It is not normally called directly by the user.

  A signed measure is a classical mathematical object
  (Diestel and Uhl, 1977)
  which can be visualised as a collection of electric charges, positive and/or
  negative, spread over the plane. Electric charges may be
  concentrated at specific points (atoms), or spread diffusely over a
  region. 

  An object of class \code{"msr"} represents a signed (i.e. real-valued)
  or vector-valued measure in the \pkg{spatstat} package.

  Spatial residuals for point process models
  (Baddeley et al, 2005, 2008) take the form of a real-valued
  or vector-valued measure. The function
  \code{\link{residuals.ppm}} returns an object of
  class \code{"msr"} representing the residual measure.

  The function \code{msr}  would not normally be called directly by the
  user. It is the low-level creator function that
  makes an object of class \code{"msr"} from raw data.
  
  The first argument \code{qscheme} is a quadrature scheme (object of
  class \code{"quad"}). It is typically created by \code{\link{quadscheme}} or
  extracted from a fitted point process model using
  \code{\link{quad.ppm}}. A quadrature scheme contains both data points
  and dummy points. The data points of \code{qscheme} are used as the locations
  of the atoms of the measure. All quadrature points
  (i.e. both data points and dummy points)
  of \code{qscheme} are used as sampling points for the density
  of the continuous component of the measure.

  The argument \code{discrete} gives the values of the
  atomic component of the measure for each \emph{data point} in \code{qscheme}.
  It should be either a numeric vector with one entry for each
  data point, or a numeric matrix with one row
  for each data point. 

  The argument \code{density} gives the values of the \emph{density}
  of the diffuse component of the measure, at each
  \emph{quadrature point} in \code{qscheme}.
  It should be either a numeric vector with one entry for each
  quadrature point, or a numeric matrix with one row
  for each quadrature point. 

  If both \code{discrete} and \code{density} are vectors
  (or one-column matrices) then the result is a signed (real-valued) measure.
  Otherwise, the result is a vector-valued measure, with the dimension
  of the vector space being determined by the number of columns
  in the matrices \code{discrete} and/or \code{density}.
  (If one of these is a \eqn{k}-column matrix and the other
  is a 1-column matrix, then the latter is replicated to \eqn{k} columns).
  
  The class \code{"msr"} has methods for \code{print},
  \code{plot} and \code{[}. 
  There is also a function \code{\link{Smooth.msr}} for smoothing a measure.
}
\value{
  An object of class \code{"msr"} that can be plotted
  by \code{\link{plot.msr}}.
}
\references{
  Baddeley, A., Turner, R., \ifelse{latex}{\out{M\o ller}}{Moller}, J.
  and Hazelton, M. (2005)
  Residual analysis for spatial point processes.
  \emph{Journal of the Royal Statistical Society, Series B}
  \bold{67}, 617--666.

  Baddeley, A., \ifelse{latex}{\out{M\o ller}}{Moller}, J.
  and Pakes, A.G. (2008) 
  Properties of residuals for spatial point processes.
  \emph{Annals of the Institute of Statistical Mathematics}
  \bold{60}, 627--649.
  
  Diestel, J. and Uhl, J.J. Jr (1977)
  \emph{Vector measures}.
  Providence, RI, USA: American Mathematical Society.
}
\author{
  \adrian
}
\seealso{
  \code{\link{plot.msr}},
  \code{\link{Smooth.msr}},
  \code{\link{[.msr}},
  \code{\link{with.msr}},
  \code{\link{split.msr}}
}
\examples{
   X <- rpoispp(function(x,y) { exp(3+3*x) })
   fit <- ppm(X, ~x+y)
   
   rp <- residuals(fit, type="pearson")
   rp

   rs <- residuals(fit, type="score")
   rs
   colnames(rs)

   # An equivalent way to construct the Pearson residual measure by hand
   Q <- quad.ppm(fit)
   lambda <- fitted(fit)
   slam <- sqrt(lambda)
   Z <- is.data(Q)
   m <- msr(Q, discrete=1/slam[Z], density = -slam)
   m
}
\keyword{spatial}
\keyword{models}
