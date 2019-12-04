\name{Gres}
\Rdversion{1.1}
\alias{Gres}
\title{
  Residual G Function
}
\description{
  Given a point process model fitted to a point pattern dataset,
  this function computes the residual \eqn{G} function,
  which serves as a diagnostic for goodness-of-fit of the model.
}
\usage{
   Gres(object, ...)
}
\arguments{
  \item{object}{
    Object to be analysed.
    Either a fitted point process model (object of class \code{"ppm"}),
    a point pattern (object of class \code{"ppp"}),
    a quadrature scheme (object of class \code{"quad"}),
    or the value returned by a previous call to \code{\link{Gcom}}.
  }
  \item{\dots}{
    Arguments passed to \code{\link{Gcom}}.
  }
}
\details{
  This command provides a diagnostic for the goodness-of-fit of
  a point process model fitted to a point pattern dataset.
  It computes a residual version of the \eqn{G} function of the
  dataset, which should be approximately zero if the model is a good
  fit to the data.

  In normal use, \code{object} is a fitted point process model
  or a point pattern. Then \code{Gres} first calls \code{\link{Gcom}}
  to compute both the nonparametric estimate of the \eqn{G} function
  and its model compensator. Then \code{Gres} computes the
  difference between them, which is the residual \eqn{G}-function.
  
  Alternatively, \code{object} may be a function value table
  (object of class \code{"fv"}) that was returned by
  a previous call to \code{\link{Gcom}}. Then \code{Gres} computes the
  residual from this object. 
}
\value{
  A function value table (object of class \code{"fv"}),
  essentially a data frame of function values.
  There is a plot method for this class. See \code{\link{fv.object}}.
}
\references{
  Baddeley, A., Rubak, E. and \Moller, J. (2011)
  Score, pseudo-score and residual
  diagnostics for spatial point process models.
  \emph{Statistical Science} \bold{26}, 613--646.
}
\author{
  \adrian,
  \ege and Jesper \Moller.
}
\seealso{
  Related functions:
  \code{\link{Gcom}}, 
  \code{\link{Gest}}.

  Alternative functions:
  \code{\link{Kres}}, 
  \code{\link{psstA}}, 
  \code{\link{psstG}}, 
  \code{\link{psst}}.

  Model-fitting:
  \code{\link{ppm}}.
}
\examples{
    data(cells)
    fit0 <- ppm(cells, ~1) # uniform Poisson
    G0 <- Gres(fit0)
    plot(G0)
# Hanisch correction estimate
    plot(G0, hres ~ r)
# uniform Poisson is clearly not correct

    fit1 <- ppm(cells, ~1, Strauss(0.08))
    plot(Gres(fit1), hres ~ r)
# fit looks approximately OK; try adjusting interaction distance

    plot(Gres(cells, interaction=Strauss(0.12)))

# How to make envelopes
    \dontrun{
    E <- envelope(fit1, Gres, model=fit1, nsim=39)
    plot(E)
    }
# For computational efficiency
    Gc <- Gcom(fit1)
    G1 <- Gres(Gc)
}
\keyword{spatial}
\keyword{models}

