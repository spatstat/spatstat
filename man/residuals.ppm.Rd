\name{residuals.ppm}
\alias{residuals.ppm}
\title{
  Residuals for Fitted Point Process Model
}
\description{
  Given a point process model fitted to a point pattern,
  compute residuals.
}
\usage{
  \method{residuals}{ppm}(object, type="raw", \dots, check=TRUE, drop=FALSE,
                fittedvalues=fitted.ppm(object, check=check, drop=drop),
                new.coef=NULL, quad=NULL)
}
\arguments{
  \item{object}{
    The fitted point process model (an object of class \code{"ppm"})
    for which residuals should be calculated.
  }
  \item{type}{
    String indicating the type of residuals to be calculated.
    Current options are
    \code{"raw"}, \code{"inverse"}, \code{"pearson"} and \code{"score"}.
    A partial match is adequate.
  }
  \item{\dots}{
    Ignored.
  }
  \item{check}{
    Logical value indicating whether to check the internal format
    of \code{object}. If there is any possibility that this object
    has been restored from a dump file, or has otherwise lost track of
    the environment where it was originally computed, set
    \code{check=TRUE}. 
  }
  \item{drop}{
    Logical value determining whether to delete quadrature points
    that were not used to fit the model. See \code{\link{quad.ppm}} for
    explanation. 
  }
  \item{fittedvalues}{
    Vector of fitted values for the conditional intensity at 
    the quadrature points,
    from which the residuals will be computed. 
    For expert use only.
  }
  \item{new.coef}{
    Optional. Numeric vector of coefficients for the model,
    replacing \code{coef(object)}.
    See the section on Modified Residuals below.
  }
  \item{quad}{
    Optional. Data specifying how to re-fit the model.
    A list of arguments passed to \code{\link{quadscheme}}.
    See the section on Modified Residuals below.
  }
}
\value{
  An object of class \code{"msr"} 
  representing a signed measure or vector-valued measure
  (see \code{\link{msr}}). This object can be plotted.
}
\details{
  This function computes several kinds of residuals for the fit of
  a point process model to a spatial point pattern dataset
  (Baddeley et al, 2005).
  Use \code{\link{plot.msr}} to plot the residuals directly,
  or \code{\link{diagnose.ppm}}
  to produce diagnostic plots based on these residuals.

  The argument \code{object} must be a fitted point process model
  (object of class \code{"ppm"}). Such objects are produced by the maximum
  pseudolikelihood fitting algorithm \code{\link{ppm}}.
  This fitted model object contains complete
  information about the original data pattern. 

  Residuals are attached both to the data points and to some
  other points in the window of observation (namely, to the dummy
  points of the quadrature scheme used to fit the model).
  If the fitted model is correct, then the sum of the 
  residuals over all (data and dummy) points in a spatial region \eqn{B}
  has mean zero. For further explanation, see Baddeley et al (2005).
  
  The type of residual
  is chosen by the argument \code{type}. Current options are

  \describe{
    \item{\code{"raw"}:}{
      the raw residuals
      \deqn{
	r_j = z_j - w_j \lambda_j
      }{
	r[j] = z[j] - w[j] lambda[j]
      }
      at the quadrature points \eqn{u_j}{u[j]},
      where \eqn{z_j}{z[j]} is the indicator equal to 1 if \eqn{u_j}{u[j]}
      is a data point and 0 if \eqn{u_j}{u[j]} is a dummy point;
      \eqn{w_j}{w[j]} is the quadrature weight attached to
      \eqn{u_j}{u[j]}; and
      \deqn{\lambda_j = \hat\lambda(u_j,x)}{lambda[j] = lambda(u[j],x)}
      is the conditional intensity of the fitted model at \eqn{u_j}{u[j]}.
      These are the spatial analogue of the martingale residuals
      of a one-dimensional counting process.
    }
    \item{\code{"inverse"}:}{
      the `inverse-lambda' residuals (Baddeley et al, 2005)
      \deqn{
	r^{(I)}_j = \frac{r_j}{\lambda_j}
        = \frac{z_j}{\lambda_j} - w_j
      }{
	rI[j] = r[j]/lambda[j] = z[j]/lambda[j] - w[j]
      }
      obtained by dividing the raw residuals by 
      the fitted conditional intensity. These are 
      a counterpart of the exponential energy marks (see \code{\link{eem}}).
    }
    \item{\code{"pearson"}:}{
      the Pearson residuals (Baddeley et al, 2005)
      \deqn{
	r^{(P)}_j = \frac{r_j}{\sqrt{\lambda_j}}
        = \frac{z_j}{\sqrt{\lambda_j}}
        - w_j \sqrt{\lambda_j}
      }{
	rP[j] = r[j]/sqrt(lambda[j])
	= z[j]/sqrt(lambda[j]) - w[j] sqrt(lambda[j])
      }
      obtained by dividing the raw residuals by the
      square root of the fitted conditional intensity.
      The Pearson residuals are standardised, in the sense
      that if the model (true and fitted) is Poisson,
      then the sum of the Pearson residuals in a spatial region \eqn{B}
      has variance equal to the area of \eqn{B}.
    }
    \item{\code{"score"}:}{
      the score residuals (Baddeley et al, 2005)
      \deqn{
	r_j = (z_j - w_j \lambda_j) x_j
      }{
	r[j] = (z[j] - w[j] lambda[j]) * x[j,]
      }
      obtained by multiplying the raw residuals \eqn{r_j}{r[j]}
      by the covariates \eqn{x_j}{x[j,]} for quadrature point \eqn{j}.
      The score residuals always sum to zero.
    }
  }
  The result of \code{residuals.ppm} is a measure
  (object of class \code{"msr"}).
  Use \code{\link{plot.msr}} to plot the residuals directly,
  or \code{\link{diagnose.ppm}} to produce diagnostic plots
  based on these residuals.
  Use \code{\link{integral.msr}} to compute the total residual.

  By default, 
  the window of the measure is the same as the original window
  of the data. If \code{drop=TRUE} then the window is the
  domain of integration of the pseudolikelihood or composite likelihood.
  This only matters when the model \code{object} was fitted using
  the border correction: in that case, if \code{drop=TRUE} the
  window of the residuals is the erosion of the original data window
  by the border correction distance \code{rbord}.
}
\section{Modified Residuals}{
  Sometimes we want to modify the calculation of residuals by using
  different values for the model parameters. This capability is
  provided by the arguments \code{new.coef} and \code{quad}.

  If \code{new.coef} is given, then the residuals will be computed
  by taking the model parameters to be \code{new.coef}.
  This should be a numeric vector
  of the same length as the vector of fitted model parameters
  \code{coef(object)}. 

  If \code{new.coef} is missing and \code{quad} is given,
  then the model parameters will
  be determined by re-fitting the model using a new
  quadrature scheme specified by \code{quad}. 
  Residuals will be computed for the
  original model \code{object} using these new parameter values.

  The argument \code{quad} should normally be
  a list of arguments in \code{name=value} format that will be
  passed to \code{\link{quadscheme}} (together with
  the original data points) to determine the new quadrature scheme.
  It may also be a quadrature scheme (object of class
  \code{"quad"}) to which the model should be fitted, or a
  point pattern (object of class \code{"ppp"}) specifying the
  \emph{dummy points} in a new quadrature scheme.
}
\references{
  Baddeley, A., Turner, R.,
  \ifelse{latex}{\out{M\o ller}}{Moller}, J. and Hazelton, M. (2005)
  Residual analysis for spatial point processes.
  \emph{Journal of the Royal Statistical Society, Series B}
  \bold{67}, 617--666.

  Baddeley, A.,
  \ifelse{latex}{\out{M\o ller}}{Moller}, J. and Pakes, A.G. (2008) 
  Properties of residuals for spatial point processes.
  \emph{Annals of the Institute of Statistical Mathematics}
  \bold{60}, 627--649.
}
\seealso{
 \code{\link{msr}},
 \code{\link{diagnose.ppm}},
 \code{\link{ppm.object}},
 \code{\link{ppm}}
}
\examples{
   fit <- ppm(cells, ~x, Strauss(r=0.15))

   # Pearson residuals
   rp <- residuals(fit, type="pe")
   rp

   # simulated data
   X <- rStrauss(100,0.7,0.05)
   # fit Strauss model 
   fit <- ppm(X, ~1, Strauss(0.05))
   res.fit <- residuals(fit)

   # check that total residual is 0 
   integral.msr(residuals(fit, drop=TRUE))

   # true model parameters
   truecoef <- c(log(100), log(0.7))
   res.true <- residuals(fit, new.coef=truecoef)  
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{models}
\keyword{methods}
