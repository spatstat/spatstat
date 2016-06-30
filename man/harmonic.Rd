\name{harmonic}
\alias{harmonic}
\title{Basis for Harmonic Functions}
\description{
  Evaluates a basis for the harmonic polynomials in \eqn{x} and \eqn{y}
  of degree less than or equal to \eqn{n}.
}
\usage{
   harmonic(x, y, n)
}
\arguments{
  \item{x}{
    Vector of \eqn{x} coordinates
  }
  \item{y}{
    Vector of \eqn{y} coordinates
  }
  \item{n}{
    Maximum degree of polynomial
  }
}
\value{
  A data frame with \code{2 * n} columns giving the values of the
  basis functions at the coordinates. Each column is labelled by an
  algebraic expression for the corresponding basis function.
}
\details{
  This function computes a basis for the harmonic polynomials
  in two variables \eqn{x} and \eqn{y} up to a given degree \eqn{n}
  and evaluates them at given \eqn{x,y} locations.
  It can be used in model formulas (for example in
  the model-fitting functions
  \code{\link{lm},\link{glm},\link{gam}} and \code{\link{ppm}}) to specify a
  linear predictor which is a harmonic function.

  A function \eqn{f(x,y)} is harmonic if
  \deqn{\frac{\partial^2}{\partial x^2} f
    + \frac{\partial^2}{\partial y^2}f = 0.}{
    (d/dx)^2 f + (d/dy)^2 f = 0.}
  The harmonic polynomials of degree less than or equal to
  \eqn{n} have a basis consisting of \eqn{2 n} functions.

  This function was implemented on a suggestion of P. McCullagh
  for fitting nonstationary spatial trend to point process models.
}
\seealso{
  \code{\link{ppm}}
}
\examples{
   # inhomogeneous point pattern
   X <- unmark(longleaf)
   \testonly{
     # smaller dataset
     X <- X[seq(1,npoints(X), by=50)]
   }

   # fit Poisson point process with log-cubic intensity
   fit.3 <- ppm(X ~ polynom(x,y,3), Poisson())

   # fit Poisson process with log-cubic-harmonic intensity
   fit.h <- ppm(X ~ harmonic(x,y,3), Poisson())

   # Likelihood ratio test
   lrts <- 2 * (logLik(fit.3) - logLik(fit.h))
   df <- with(coords(X),
              ncol(polynom(x,y,3)) - ncol(harmonic(x,y,3)))
   pval <- 1 - pchisq(lrts, df=df)
}
\author{
  \spatstatAuthors.
}
\keyword{spatial}
\keyword{models}
