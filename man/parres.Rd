\name{parres}
\alias{parres}
\title{
  Partial Residuals for Point Process Model
}
\description{
  Computes the smoothed partial residuals, a diagnostic
  for transformation of a covariate in a Poisson point process model.
}
\usage{
parres(model, covariate, ...,
       smooth.effect=FALSE, subregion=NULL,
       bw = "nrd0", adjust=1, from = NULL, to = NULL, n = 512,
       bw.input = c("points", "quad"), bw.restrict=FALSE, covname)
}
\arguments{
  \item{model}{
    Fitted point process model (object of class \code{"ppm"}).
  }
  \item{covariate}{
    The covariate of interest.
    Either a character string matching the name of one of
    the canonical covariates in the model, or one of the names
    \code{"x"} or \code{"y"} referring to the Cartesian coordinates,
    or one of the names of the covariates given when \code{model} was fitted,
    or a pixel image (object of class \code{"im"})
    or \code{function(x,y)} supplying the
    values of a covariate at any location.
  }
  \item{smooth.effect}{
    Logical. Determines the choice of algorithm. See Details.
  }
  \item{subregion}{
    Optional.  A window (object of class \code{"owin"})
    specifying a subset of the spatial domain of the data.
    The calculation will be confined to the data in this subregion.
  }
  \item{bw}{
    Smoothing bandwidth or bandwidth rule
    (passed to \code{\link[stats]{density.default}}).
  }
  \item{adjust}{
    Smoothing bandwidth adjustment factor
    (passed to \code{\link[stats]{density.default}}).
  }
  \item{n, from, to}{
    Arguments passed to \code{\link[stats]{density.default}} to
    control the number and range of values at which the function
    will be estimated.
  }
  \item{\dots}{
    Additional arguments passed to \code{\link[stats]{density.default}}.
  }
  \item{bw.input}{
    Character string specifying the input data used for automatic
    bandwidth selection. 
  }
  \item{bw.restrict}{
    Logical value, specifying whether bandwidth selection is performed using
    data from the entire spatial domain or from the \code{subregion}.
  }
  \item{covname}{
    Optional. Character string to use as the name of the covariate.
  }
}
\details{
  This command computes the smoothed partial residual diagnostic
  (Baddeley, Chang, Song and Turner, 2012) 
  for the transformation of a covariate
  in a Poisson point process model.
  
  The argument \code{model} must be a fitted Poisson point process model.
  
  The diagnostic works in two different ways:
  \describe{
    \item{Canonical covariate:}{
      The argument \code{covariate} may be a character string
      which is the name of one of the \emph{canonical covariates} in the
      model. 
      The canonical covariates are the
      functions \eqn{Z_j}{Z[j]} that appear
      in the expression for the Poisson point process intensity
      \deqn{
	\lambda(u) = \exp(\beta_1 Z_1(u) + \ldots + \beta_p Z_p(u))
      }{
	lambda(u) = exp(beta[1] * Z[1](u) + \ldots + \beta[p] * Z[p](u))
      }
      at spatial location \eqn{u}.
      Type \code{names(coef(model))} to see the names of the
      canonical covariates in \code{model}.
      If the selected covariate is \eqn{Z_j}{Z[j]}, then 
      the diagnostic plot concerns the model term
      \eqn{\beta_j Z_j(u)}{beta[j] * Z[j](u)}. The plot shows a smooth
      estimate of a function \eqn{h(z)} that should replace this linear
      term, that is, \eqn{\beta_j Z_j(u)}{beta[j] * Z[j](u)} should be
      replaced by \eqn{h(Z_j(u))}{h(Z[j](u))}. The linear function is
      also plotted as a dotted line.
    }
    \item{New covariate:}{
      If the argument \code{covariate} is a pixel image
      (object of class \code{"im"}) or a \code{function(x,y)},
      it is assumed to provide the values of a covariate that is
      not present in the model.
      Alternatively \code{covariate} can be the name of a
      covariate that was supplied when the model was fitted
      (i.e. in the call to \code{\link{ppm}})
      but which does not feature in the model formula.
      In either case we speak of a new covariate \eqn{Z(u)}.
      If the fitted model intensity is \eqn{\lambda(u)}{lambda(u)}
      then we consider modifying this to
      \eqn{\lambda(u) \exp(h(Z(u)))}{lambda(u) * exp(h(Z(u)))}
      where \eqn{h(z)} is some function. The diagnostic plot shows
      an estimate of \eqn{h(z)}.
      \bold{Warning: in this case the diagnostic is not theoretically
      justified. This option is provided for research purposes.}
    }
  }
  Alternatively \code{covariate} can be one of the character strings
  \code{"x"} or \code{"y"} signifying the Cartesian coordinates.
  The behaviour here depends on whether the coordinate was one of the
  canonical covariates in the model.

  If there is more than one canonical covariate in the model
  that depends on the specified \code{covariate}, then
  the covariate effect is computed using all these canonical covariates.
  For example in a log-quadratic model which includes the terms \code{x} and
  \code{I(x^2)}, the quadratic effect involving both these terms
  will be computed.

  There are two choices for the algorithm.
  If \code{smooth.effect=TRUE}, the fitted covariate effect (according
  to \code{model}) is added to the point process residuals, then
  smoothing is applied to these values. If \code{smooth.effect=FALSE},
  the point process residuals are smoothed first, and then the fitted
  covariate effect is added to the result.

  The smoothing bandwidth is controlled by the arguments
  \code{bw}, \code{adjust}, \code{bw.input} and \code{bw.restrict}.
  If \code{bw} is a numeric value, then
  the bandwidth is taken to be \code{adjust * bw}.
  If \code{bw} is a string representing a bandwidth selection rule
  (recognised by \code{\link[stats]{density.default}})
  then the bandwidth is selected by this rule.

  The data used for automatic bandwidth selection are
  specified by \code{bw.input} and \code{bw.restrict}.
  If \code{bw.input="points"}  (the default) then bandwidth selection is
  based on the covariate values at the points of the original point
  pattern dataset to which the model was fitted.
  If \code{bw.input="quad"} then bandwidth selection is
  based on the covariate values at every quadrature point used to
  fit the model.
  If \code{bw.restrict=TRUE} then the bandwidth selection is performed
  using only data from inside the \code{subregion}.
}
\section{Slow computation}{
  In a large dataset, computation can be very slow if the default
  settings are used, because the smoothing bandwidth is selected
  automatically. To avoid this, specify a numerical value
  for the bandwidth \code{bw}. One strategy is to use a coarser
  subset of the data to select \code{bw} automatically.
  The selected bandwidth can be read off the print output for
  \code{parres}.  
}
\value{
  A function value table (object of class \code{"fv"})
  containing the values of the smoothed partial residual,
  the estimated variance, and the fitted effect of the covariate.
  Also belongs to the class \code{"parres"}
  which has methods for \code{print} and \code{plot}.
}
\references{
  Baddeley, A., Chang, Y.-M., Song, Y. and Turner, R. (2013)
  Residual diagnostics for covariate effects in
  spatial point process models.
  \emph{Journal of Computational and Graphical Statistics},
  \bold{22}, 886--905.
}
\author{
  Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/},
  Rolf Turner
  \email{r.turner@auckland.ac.nz},
  Ya-Mei Chang and Yong Song.
}
\seealso{
  \code{\link{addvar}},
  \code{\link{rhohat}},
  \code{\link{rho2hat}}
}
\examples{
  X <-  rpoispp(function(x,y){exp(3+x+2*x^2)})
  model <- ppm(X, ~x+y)
  tra <- parres(model, "x")
  plot(tra)
  plot(parres(model, "x", subregion=square(0.5)))
  model2 <- ppm(X, ~x+I(x^2)+y)
  plot(parres(model2, "x"))
  Z <- setcov(owin())
  plot(parres(model2, Z))
}
\keyword{spatial}
\keyword{models}
