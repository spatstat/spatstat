\name{Kres}
\Rdversion{1.1}
\alias{Kres}
\title{
  Residual K Function
}
\description{
  Given a point process model fitted to a point pattern dataset,
  this function computes the residual \eqn{K} function,
  which serves as a diagnostic for goodness-of-fit of the model.
}
\usage{
   Kres(object, ...)
}
\arguments{
  \item{object}{
    Object to be analysed.
    Either a fitted point process model (object of class \code{"ppm"}),
    a point pattern (object of class \code{"ppp"}),
    a quadrature scheme (object of class \code{"quad"}),
    or the value returned by a previous call to \code{\link{Kcom}}.
  }
  \item{\dots}{
    Arguments passed to \code{\link{Kcom}}.
  }
}
\details{
  This command provides a diagnostic for the goodness-of-fit of
  a point process model fitted to a point pattern dataset.
  It computes a residual version of the \eqn{K} function of the
  dataset, which should be approximately zero if the model is a good
  fit to the data.

  In normal use, \code{object} is a fitted point process model
  or a point pattern. Then \code{Kres} first calls \code{\link{Kcom}}
  to compute both the nonparametric estimate of the \eqn{K} function
  and its model compensator. Then \code{Kres} computes the
  difference between them, which is the residual \eqn{K}-function.
  
  Alternatively, \code{object} may be a function value table
  (object of class \code{"fv"}) that was returned by
  a previous call to \code{\link{Kcom}}. Then \code{Kres} computes the
  residual from this object. 
}
\value{
  A function value table (object of class \code{"fv"}),
  essentially a data frame of function values.
  There is a plot method for this class. See \code{\link{fv.object}}.
}
\references{
  Baddeley, A., Rubak, E. and Moller, J. (2011)
  Score, pseudo-score and residual
  diagnostics for spatial point process models.
  \emph{Statistical Science} \bold{26}, 613--646.
}
\author{
  Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  Ege Rubak and Jesper Moller.
}
\seealso{
  Related functions:
  \code{\link{Kcom}}, 
  \code{\link{Kest}}.
  
  Alternative functions:
  \code{\link{Gres}},
  \code{\link{psstG}}, \code{\link{psstA}},  \code{\link{psst}}.

  Point process models: \code{\link{ppm}}.
}
\examples{
    data(cells)
    fit0 <- ppm(cells, ~1) # uniform Poisson
    \testonly{    fit0 <- ppm(cells, ~1, nd=16)}
    K0 <- Kres(fit0)
    K0
    plot(K0)
# isotropic-correction estimate
    plot(K0, ires ~ r)
# uniform Poisson is clearly not correct

    fit1 <- ppm(cells, ~1, Strauss(0.08))
    \testonly{fit1 <- ppm(cells, ~1, Strauss(0.08), nd=16)}
    K1 <- Kres(fit1)

    if(interactive()) {
      plot(K1, ires ~ r)
   # fit looks approximately OK; try adjusting interaction distance
      plot(Kres(cells, interaction=Strauss(0.12)))
    }

# How to make envelopes
    \dontrun{
    E <- envelope(fit1, Kres, interaction=as.interact(fit1), nsim=19)
    plot(E)
    }

# For computational efficiency
    Kc <- Kcom(fit1)
    K1 <- Kres(Kc)
}
\keyword{spatial}
\keyword{models}

