\name{psst}
\alias{psst}
\title{
  Pseudoscore Diagnostic For Fitted Model against General Alternative
}
\description{
  Given a point process model fitted to a point pattern dataset,
  and any choice of functional summary statistic,
  this function computes the pseudoscore test statistic
  of goodness-of-fit for the model.
}
\usage{
psst(object, fun, r = NULL, breaks = NULL, ...,
     model=NULL,
     trend = ~1, interaction = Poisson(), rbord = reach(interaction),
     truecoef=NULL, hi.res=NULL, funargs = list(correction="best"),
     verbose=TRUE)
}
\arguments{
  \item{object}{
    Object to be analysed.
    Either a fitted point process model (object of class \code{"ppm"})
    or a point pattern (object of class \code{"ppp"})
    or quadrature scheme (object of class \code{"quad"}).
  }
  \item{fun}{
    Summary function to be applied to each point pattern.
  }
  \item{r}{
    Optional. 
    Vector of values of the argument \eqn{r} at which the
    function \eqn{S(r)} should be computed.
    This argument is usually not specified. There is a sensible default.
  }
  \item{breaks}{
    Optional alternative to \code{r} for advanced use.
  }
  \item{\dots}{
    Ignored.
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
  \item{funargs}{
    List of additional arguments to be passed to \code{fun}.
  }
  \item{verbose}{
    Logical value determining whether to print progress reports
    during the computation.
  }
}
\details{
  Let \eqn{x} be a point pattern dataset consisting of points
  \eqn{x_1,\ldots,x_n}{x[1],...,x[n]} in a window \eqn{W}.
  Consider a point process model fitted to \eqn{x}, with
  conditional intensity
  \eqn{\lambda(u,x)}{lambda(u,x)} at location \eqn{u}.
  For the purpose of testing goodness-of-fit, we regard the fitted model
  as the null hypothesis. Given a functional summary statistic \eqn{S},
  consider a family of alternative models obtained by exponential
  tilting of the null model by \eqn{S}. 
  The pseudoscore for the null model is
  \deqn{
    V(r) = \sum_i \Delta S(x_i, x, r) - \int_W \Delta S(u,x, r) \lambda(u,x)
    {\rm d} u
  }{
    V(r) = sum( Delta S(x[i], x, r)) - integral( Delta S(u,x, r) lambda(u,x) du)
  }
  where the \eqn{\Delta}{Delta} operator is
  \deqn{
    \Delta S(u,x, r) = S(x\cup\{u\}, r) - S(x\setminus u, r)
  }{
    Delta S(u,x, r) = S(x union u, r) - S(x setminus u, r)
  }
  the difference between the values of \eqn{S} for the
  point pattern with and without the point \eqn{u}.

  According to the Georgii-Nguyen-Zessin formula, \eqn{V(r)} should have
  mean zero if the model is correct (ignoring the fact that the
  parameters of the model have been estimated). Hence \eqn{V(r)} can be
  used as a diagnostic for goodness-of-fit.

  This algorithm computes \eqn{V(r)} by direct evaluation of the sum and
  integral. It is computationally intensive, but it is available for
  any summary statistic \eqn{S(r)}.

  The diagnostic \eqn{V(r)} is also called 
  the \bold{pseudoresidual} of \eqn{S}. On the right
  hand side of the equation for \eqn{V(r)} given above,
  the sum over points of \eqn{x} is called the
  \bold{pseudosum} and the integral is called the \bold{pseudocompensator}.
}
\value{
  A function value table (object of class \code{"fv"}),
  essentially a data frame of function values.

  Columns in this data frame include \code{dat} for the pseudosum,
  \code{com} for the compensator and \code{res} for the
  pseudoresidual.
  
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
  Special cases:
  \code{\link{psstA}},
  \code{\link{psstG}}.

  Alternative functions:
  \code{\link{Kres}},
  \code{\link{Gres}}.
}
\examples{
    data(cells)
    fit0 <- ppm(cells, ~1) # uniform Poisson
    \testonly{fit0 <- ppm(cells, ~1, nd=8)}
    G0 <- psst(fit0, Gest)
    G0
    if(interactive()) plot(G0)
}
\keyword{spatial}
\keyword{models}
