\name{Kcom}
\Rdversion{1.1}
\alias{Kcom}
\title{
  Model Compensator of K Function
}
\description{
  Given a point process model fitted to a point pattern dataset,
  this function computes the \emph{compensator} 
  of the \eqn{K} function based on the fitted model 
  (as well as the usual nonparametric estimates
  of \eqn{K} based on the data alone).
  Comparison between the nonparametric and model-compensated \eqn{K}
  functions serves as a diagnostic for the model.
}
\usage{
Kcom(object, r = NULL, breaks = NULL, ...,
     correction = c("border", "isotropic", "translate"),
     conditional = !is.poisson(object),
     restrict = FALSE,
     model = NULL,
     trend = ~1, interaction = Poisson(), rbord = reach(interaction),
     compute.var = TRUE,
     truecoef = NULL, hi.res = NULL)
}
\arguments{
  \item{object}{
    Object to be analysed.
    Either a fitted point process model (object of class \code{"ppm"})
    or a point pattern (object of class \code{"ppp"})
    or quadrature scheme (object of class \code{"quad"}).
}
  \item{r}{
    Optional. 
    Vector of values of the argument \eqn{r} at which the
    function \eqn{K(r)} should be computed.
    This argument is usually not specified. There is a sensible default.
  }
  \item{breaks}{
	This argument is for advanced use only.
  }
  \item{\dots}{
    Ignored.
  }
  \item{correction}{
    Optional vector of character strings specifying the edge
    correction(s) to be used. See \code{\link{Kest}} for options.
  }
  \item{conditional}{
    Optional. Logical value indicating whether to 
    compute the estimates for the conditional case. See Details.
  }
  \item{restrict}{
    Logical value indicating
    whether to compute the restriction estimator (\code{restrict=TRUE}) or
    the reweighting estimator (\code{restrict=FALSE}, the default).
    Applies only if \code{conditional=TRUE}.
    See Details.
  }
  \item{model}{
    Optional. A fitted point process model (object of
    class \code{"ppm"}) to be re-fitted to the data
    using \code{\link{update.ppm}}, if \code{object} is a point pattern.
    Overrides the arguments \code{trend,interaction,rbord}.
  }
  \item{trend,interaction,rbord}{
    Optional. Arguments passed to \code{\link{ppm}}
    to fit a point process model to the data,
    if \code{object} is a point pattern.
    See \code{\link{ppm}} for details.
  }
  \item{compute.var}{
    Logical value indicating whether to compute the
    Poincare variance bound for the residual \eqn{K} function
    (calculation is only implemented for the isotropic correction).
  }
  \item{truecoef}{
    Optional. Numeric vector. If present, this will be treated as 
    if it were the true coefficient vector of the point process model,
    in calculating the diagnostic. Incompatible with \code{hi.res}.
  }
  \item{hi.res}{
    Optional. List of parameters passed to \code{\link{quadscheme}}.
    If this argument is present, the model will be
    re-fitted at high resolution as specified by these parameters.
    The coefficients
    of the resulting fitted model will be taken as the true coefficients.
    Then the diagnostic will be computed for the default
    quadrature scheme, but using the high resolution coefficients.
  }
}
\details{
  This command provides a diagnostic for the goodness-of-fit of
  a point process model fitted to a point pattern dataset.
  It computes an estimate of the \eqn{K} function of the
  dataset, together with a \emph{model compensator} of the
  \eqn{K} function, which should be approximately equal if the model is a good
  fit to the data.

  The first argument, \code{object}, is usually a fitted point process model
  (object of class \code{"ppm"}), obtained from the
  model-fitting function \code{\link{ppm}}.

  For convenience, \code{object} can also be a point pattern
  (object of class \code{"ppp"}). In that case, a point process
  model will be fitted to it, by calling \code{\link{ppm}} using the arguments
  \code{trend} (for the first order trend),
  \code{interaction} (for the interpoint interaction)
  and \code{rbord} (for the erosion distance in the border correction
  for the pseudolikelihood). See \code{\link{ppm}} for details
  of these arguments.

  The algorithm first extracts the original point pattern dataset
  (to which the model was fitted) and computes the 
  standard nonparametric estimates of the \eqn{K} function.
  It then also computes the \emph{model compensator} of the
  \eqn{K} function. The different function estimates are returned
  as columns in a data frame (of class \code{"fv"}).

  The argument \code{correction} determines the edge correction(s)
  to be applied. See \code{\link{Kest}} for explanation of the principle
  of edge corrections. The following table gives the options
  for the \code{correction} argument, and the corresponding
  column names in the result:

  \tabular{llll}{
    \code{correction} \tab
    \bold{description of correction} \tab
    \bold{nonparametric} \tab \bold{compensator} \cr
    \code{"isotropic"} \tab Ripley isotropic correction
    \tab \code{iso} \tab \code{icom} \cr
    \code{"translate"} \tab Ohser-Stoyan translation correction
    \tab \code{trans} \tab \code{tcom} \cr
    \code{"border"} \tab border correction
    \tab \code{border} \tab \code{bcom}
  }

  The nonparametric estimates can all be expressed in the form
  \deqn{
    \hat K(r) = \sum_i \sum_{j < i} e(x_i,x_j,r,x) I\{ d(x_i,x_j) \le r \}
  }{
    K(r) = sum[i] sum[j < i] e(x[i], x[j], r, x) I( d(x[i],x[j]) <= r )
  }
  where \eqn{x_i}{x[i]} is the \eqn{i}-th data point,
  \eqn{d(x_i,x_j)}{d(x[i],x[j])} is the distance between \eqn{x_i}{x[i]} and
  \eqn{x_j}{x[j]}, and \eqn{e(x_i,x_j,r,x)}{e(x[i],x[j],r,x)} is
  a term that serves to correct edge effects and to re-normalise the
  sum. The corresponding model compensator is
  \deqn{
    {\bf C} \, \tilde K(r) = \int_W \lambda(u,x) \sum_j e(u,x_j,r,x \cup u) I\{ d(u,x_j) \le r\}
  }{
    C K(r) = integral[u] lambda(u,x) sum[j] e(u, x[j], r, x+u) I( d(u,x[j]) <= r )
  }
  where the integral is over all locations \eqn{u} in
  the observation window, 
  \eqn{\lambda(u,x)}{lambda(u,x)} denotes the conditional intensity
  of the model at the location \eqn{u}, and \eqn{x \cup u}{x+u} denotes the
  data point pattern \eqn{x} augmented by adding the extra point \eqn{u}.
  
  If the fitted model is a Poisson point process, then the formulae above
  are exactly what is computed. If the fitted model is not Poisson, the 
  formulae above are modified slightly to handle edge effects.

  The modification is determined by the arguments
  \code{conditional} and \code{restrict}.
  The value of \code{conditional} defaults to \code{FALSE} for Poisson models
  and \code{TRUE} for non-Poisson models.
  If \code{conditional=FALSE} then the formulae above are not modified.
  If \code{conditional=TRUE}, then the algorithm calculates
  the \emph{restriction estimator} if \code{restrict=TRUE},
  and calculates the \emph{reweighting estimator} if \code{restrict=FALSE}.
  See Appendix D of Baddeley, Rubak and \ifelse{latex}{\out{M\o ller}}{Moller} (2011).
  Thus, by default, the reweighting estimator is computed
  for non-Poisson models.
  
  The nonparametric estimates of \eqn{K(r)} are approximately unbiased
  estimates of the \eqn{K}-function, assuming the point process is
  stationary. The model compensators are unbiased estimates
  \emph{of the mean values of the corresponding nonparametric estimates},
  assuming the model is true. Thus, if the model is a good fit, the mean value
  of the difference between the nonparametric estimates and model compensators
  is approximately zero.
}
\value{
  A function value table (object of class \code{"fv"}),
  essentially a data frame of function values.
  There is a plot method for this class. See \code{\link{fv.object}}.
}
\references{
  Baddeley, A., Rubak, E. and \ifelse{latex}{\out{M\o ller}}{Moller}, J. (2011)
  Score, pseudo-score and residual
  diagnostics for spatial point process models.
  \emph{Statistical Science} \bold{26}, 613--646.
}
\author{
  \adrian
  
  
  \ege and Jesper \ifelse{latex}{\out{M\o ller}}{Moller}.
}
\seealso{
  Related functions: 
  \code{\link{Kres}},
  \code{\link{Kest}}.

  Alternative functions:
  \code{\link{Gcom}},
  \code{\link{psstG}}, \code{\link{psstA}},  \code{\link{psst}}.

  Point process models: \code{\link{ppm}}.
}
\examples{
    fit0 <- ppm(cells, ~1) # uniform Poisson
    \testonly{fit0 <- ppm(cells, ~1, nd=16)}

    if(interactive()) {
      plot(Kcom(fit0))
# compare the isotropic-correction estimates
      plot(Kcom(fit0), cbind(iso, icom) ~ r)
# uniform Poisson is clearly not correct
    }

    fit1 <- ppm(cells, ~1, Strauss(0.08))
    \testonly{fit1 <- ppm(cells, ~1, Strauss(0.08), nd=16)}
    K1 <- Kcom(fit1)
    K1
    if(interactive()) {
      plot(K1)
      plot(K1, cbind(iso, icom) ~ r)
      plot(K1, cbind(trans, tcom) ~ r)
# how to plot the difference between nonparametric estimates and compensators
      plot(K1, iso - icom ~ r)
# fit looks approximately OK; try adjusting interaction distance
    }
    fit2 <- ppm(cells, ~1, Strauss(0.12))
    \testonly{fit2 <- ppm(cells, ~1, Strauss(0.12), nd=16)}
    K2 <- Kcom(fit2)
    if(interactive()) {
      plot(K2)
      plot(K2, cbind(iso, icom) ~ r)
      plot(K2, iso - icom ~ r)
    }
}
\keyword{spatial}
\keyword{models}
