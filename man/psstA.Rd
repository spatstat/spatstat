\name{psstA}
\Rdversion{1.1}
\alias{psstA}
\title{
  Pseudoscore Diagnostic For Fitted Model against Area-Interaction Alternative
}
\description{
  Given a point process model fitted to a point pattern dataset,
  this function computes the pseudoscore diagnostic 
  of goodness-of-fit for the model, against moderately
  clustered or moderately inhibited alternatives of area-interaction type.
}
\usage{
psstA(object, r = NULL, breaks = NULL, \dots,
      model = NULL, 
      trend = ~1, interaction = Poisson(),
      rbord = reach(interaction), ppmcorrection = "border",
      correction = "all",
      truecoef = NULL, hi.res = NULL,
      nr=spatstat.options("psstA.nr"),
      ngrid=spatstat.options("psstA.ngrid"))
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
    diagnostic should be computed.
    This argument is usually not specified. There is a sensible default.
  }
  \item{breaks}{
	This argument is for internal use only.
  }
  \item{\dots}{
    Extra arguments passed to \code{\link{quadscheme}} to determine
    the quadrature scheme, if \code{object} is a point pattern.
  }
  \item{model}{
    Optional. A fitted point process model (object of
    class \code{"ppm"}) to be re-fitted to the data
    using \code{\link{update.ppm}}, if \code{object} is a point pattern.
    Overrides the arguments \code{trend,interaction,rbord,ppmcorrection}.
  }
  \item{trend,interaction,rbord}{
    Optional. Arguments passed to \code{\link{ppm}}
    to fit a point process model to the data,
    if \code{object} is a point pattern.
    See \code{\link{ppm}} for details.
  }
  \item{ppmcorrection}{
    Optional. Character string specifying the edge correction
    for the pseudolikelihood to be used
    in fitting the point process model. Passed to \code{\link{ppm}}.
  }
  \item{correction}{
    Optional. Character string specifying which diagnostic quantities
    will be computed. Options are \code{"all"} and \code{"best"}.
    The default is to compute all diagnostic quantities.
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
  \item{nr}{
    Optional. Number of \code{r} values to be used
    if \code{r} is not specified.
  }
  \item{ngrid}{
    Integer. Number of points in the square grid used to compute
    the approximate area.
  }
}
\details{
  This function computes the pseudoscore test statistic
  which can be used as a diagnostic for goodness-of-fit of a fitted
  point process model.
  
  Let \eqn{x} be a point pattern dataset consisting of points
  \eqn{x_1,\ldots,x_n}{x[1],...,x[n]} in a window \eqn{W}.
  Consider a point process model fitted to \eqn{x}, with
  conditional intensity
  \eqn{\lambda(u,x)}{lambda(u,x)} at location \eqn{u}.
  For the purpose of testing goodness-of-fit, we regard the fitted model
  as the null hypothesis.
  The alternative hypothesis is a family of
  hybrid models obtained by combining 
  the fitted model with the area-interaction process
  (see \code{\link{AreaInter}}). The family of alternatives includes
  models that are slightly more regular than the fitted model,
  and others that are slightly more clustered than the fitted model.

  The pseudoscore, evaluated at the null model, is
  \deqn{
    V(r) = \sum_i A(x_i, x, r) - \int_W A(u,x, r) \lambda(u,x)
    {\rm d} u
  }{
    V(r) = sum( A(x[i], x, r)) - integral( A(u,x,r) lambda(u,x) du)
  }
  where
  \deqn{
    A(u,x,r) = B(x\cup\{u\},r) - B(x\setminus u, r)
  }{
    A(u,x,r) = B(x union u, r) - B(x setminus u, r)
  }
  where \eqn{B(x,r)} is the area of the union of the discs of radius
  \eqn{r} centred at the points of \eqn{x} (i.e. \eqn{B(x,r)} is the area
  of the dilation of \eqn{x} by a distance \eqn{r}).
  Thus \eqn{A(u,x,r)} is the \emph{unclaimed area} associated with
  \eqn{u}, that is, the area of that part of the disc
  of radius \eqn{r} centred at the point \eqn{u} that is
  not covered by any of the discs of radius \eqn{r} centred at
  points of \eqn{x}.

  According to the Georgii-Nguyen-Zessin formula, \eqn{V(r)} should have
  mean zero if the model is correct (ignoring the fact that the
  parameters of the model have been estimated). Hence \eqn{V(r)} can be
  used as a diagnostic for goodness-of-fit.

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
\section{Warning}{
  This computation can take a \bold{very long time}.

  To shorten the computation time, choose smaller values of the
  arguments \code{nr} and \code{ngrid}, or reduce the values of their
  defaults \code{spatstat.options("psstA.nr")}
  and \code{spatstat.options("psstA.ngrid")}.

  Computation time is roughly proportional to
  \code{nr * npoints * ngrid^2} where \code{npoints} is the number
  of points in the point pattern.
}
\references{
  Baddeley, A., Rubak, E. and \ifelse{latex}{\out{M\o ller}}{Moller}, J. (2011)
  Score, pseudo-score and residual
  diagnostics for spatial point process models.
  \emph{Statistical Science} \bold{26}, 613--646.
}
\author{
  Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  Ege Rubak and Jesper \ifelse{latex}{\out{M\o ller}}{Moller}.
}
\seealso{
  Alternative functions:
  \code{\link{psstG}},
  \code{\link{psst}},
  \code{\link{Gres}},
  \code{\link{Kres}}.

  Point process models: \code{\link{ppm}}.
  
  Options: \code{\link{spatstat.options}}
}
\examples{
   pso <- spatstat.options(psstA.ngrid=16,psstA.nr=10)
   X <- rStrauss(200,0.1,0.05)
   plot(psstA(X))
   plot(psstA(X, interaction=Strauss(0.05)))
   spatstat.options(pso)
}
\keyword{spatial}
\keyword{models}
