\name{Gcom}
\Rdversion{1.1}
\alias{Gcom}
\title{
  Model Compensator of Nearest Neighbour Function
}
\description{
  Given a point process model fitted to a point pattern dataset,
  this function computes the \emph{compensator} 
  of the nearest neighbour distance distribution function \eqn{G}
  based on the fitted model 
  (as well as the usual nonparametric estimates
  of \eqn{G} based on the data alone).
  Comparison between the nonparametric and model-compensated \eqn{G}
  functions serves as a diagnostic for the model.
}
\usage{
Gcom(object, r = NULL, breaks = NULL, ...,
     correction = c("border", "Hanisch"),
     conditional = !is.poisson(object),
     restrict=FALSE,
     model=NULL,
     trend = ~1, interaction = Poisson(),
     rbord = reach(interaction),
     ppmcorrection="border",
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
    function \eqn{G(r)} should be computed.
    This argument is usually not specified. There is a sensible default.
  }
  \item{breaks}{
    This argument is for internal use only.
  }
  \item{correction}{
    Edge correction(s) to be employed in calculating the compensator.
    Options are \code{"border"}, \code{"Hanisch"} and \code{"best"}.
    Alternatively \code{correction="all"} selects all options.
  }
  \item{conditional}{
    Optional. Logical value indicating whether to 
    compute the estimates for the conditional case. See Details.
  }
  \item{restrict}{
    Logical value indicating whether to compute
    the restriction estimator (\code{restrict=TRUE}) or
    the reweighting estimator (\code{restrict=FALSE}, the default).
    Applies only if \code{conditional=TRUE}.  See Details.
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
  \item{\dots}{
    Extra arguments passed to \code{\link{ppm}}.
  }
  \item{ppmcorrection}{
    The \code{correction} argument to \code{\link{ppm}}.
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
  It computes different estimates of the nearest neighbour distance
  distribution function \eqn{G} of the
  dataset, which should be approximately equal if the model is a good
  fit to the data.

  The first argument, \code{object}, is usually a fitted point process model
  (object of class \code{"ppm"}), obtained from the
  model-fitting function \code{\link{ppm}}.

  For convenience, \code{object} can also be a point pattern
  (object of class \code{"ppp"}).
  In that case, a point process
  model will be fitted to it,
  by calling \code{\link{ppm}} using the arguments
  \code{trend} (for the first order trend),
  \code{interaction} (for the interpoint interaction)
  and \code{rbord} (for the erosion distance in the border correction
  for the pseudolikelihood). See \code{\link{ppm}} for details
  of these arguments.
  
  The algorithm first extracts the original point pattern dataset
  (to which the model was fitted) and computes the 
  standard nonparametric estimates of the \eqn{G} function.
  It then also computes the \emph{model-compensated} 
  \eqn{G} function. The different functions are returned
  as columns in a data frame (of class \code{"fv"}).
  The interpretation of the columns is as follows
  (ignoring edge corrections):
  \describe{
    \item{\code{bord}:}{ 
      the nonparametric border-correction estimate of \eqn{G(r)},
      \deqn{
	\hat G(r) = \frac{\sum_i I\{ d_i \le r\} I\{ b_i > r \}}{\sum_i I\{
	  b_i > r\}}
      }{
	G(r) = (sum[i] I(d[i] <= r) I(b[i] > r))/(sum[i] I(b[i] > r))
      }
      where \eqn{d_i}{d[i]} is the distance from the \eqn{i}-th data point
      to its nearest neighbour, and \eqn{b_i}{b[i]} is the distance from the
      \eqn{i}-th data point to the boundary of the window \eqn{W}.
    }
    \item{\code{bcom}:}{
      the model compensator of the border-correction estimate
      \deqn{
	{\bf C}\,  \hat G(r) = \frac{\int \lambda(u,x) I\{ b(u) > r\} I\{ d(u,x)
	  \le r\}}{ 1 + \sum_i I\{ b_i > r\} }
      }{
	C G(r) = (integral[u] lambda(u,x) I(b(u) > r) I( d(u,x) <= r ))/(1
	+ sum[i] I(b[i] > r))
      }
      where 
      \eqn{\lambda(u,x)}{lambda(u,x)} denotes the conditional intensity
      of the model at the location \eqn{u}, and \eqn{d(u,x)} denotes the
      distance from \eqn{u} to the nearest point in \eqn{x}, while
      \eqn{b(u)} denotes the distance from \eqn{u} to the boundary of the
      window\eqn{W}.
    }
    \item{\code{han}:}{
      the nonparametric Hanisch estimate of \eqn{G(r)}
      \deqn{
	\hat G(r) = \frac{D(r)}{D(\infty)}
      }{
	G(r) = D(r)/D(infty)
      }
      where
      \deqn{
	D(r) = \sum_i 
	\frac{ I\{x_i \in W_{\ominus d_i}\} I\{d_i \le r\}
	}{
	  \mbox{area}(W_{\ominus d_i})
	}
      }{
	D(r) = sum[i] I(x[i] in W[-r]) I(d[i] <= r)/area(W[-d[i]])
      }
      in which \eqn{W_{\ominus r}}{W[-r]} denotes the erosion of the window
      \eqn{W} by a distance \eqn{r}.
    }
    \item{\code{hcom}:}{
      the corresponding model-compensated function 
      \deqn{
	{\bf C} \, G(r) = \int_W 
	\frac{
	  \lambda(u,x) I(u \in W_{\ominus d(u)}) I(d(u) \le r)
	}{ 
	  \hat D(\infty) \mbox{area}(W_{\ominus d(u)}) + 1
	}
      }{
	C G(r) = integral[u] lambda(u,x) I(u in W[-d(u)]) I(d(u) <= r)/
	(1 + D(infty) area(W[-d(u)]))
      }
      where \eqn{d(u) = d(u, x)} is the (`empty space') 
      distance from location \eqn{u} to the nearest point of \eqn{x}.
    }
  }
  
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
  See Appendix E of Baddeley, Rubak
  and \ifelse{latex}{\out{M\o ller}}{Moller} (2011).
  See also \code{\link{spatstat.options}('eroded.intensity')}.
  Thus, by default, the reweighting estimator is computed
  for non-Poisson models.

  The border-corrected and Hanisch-corrected estimates of \eqn{G(r)} are
  approximately unbiased estimates of the \eqn{G}-function,
  assuming the point process is
  stationary. The model-compensated functions are unbiased estimates
  \emph{of the mean value of the corresponding nonparametric estimate},
  assuming the model is true. Thus, if the model is a good fit, the mean value
  of the difference between the nonparametric and model-compensated
  estimates is approximately zero.

  To compute the difference between the nonparametric and model-compensated
  functions, use \code{\link{Gres}}.
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
  Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  Ege Rubak \email{rubak@math.aau.dk} and Jesper \ifelse{latex}{\out{M\o ller}}{Moller}.
}
\seealso{
  Related functions:
  \code{\link{Gest}},
  \code{\link{Gres}}.

  Alternative functions:
  \code{\link{Kcom}}, 
  \code{\link{psstA}}, 
  \code{\link{psstG}}, 
  \code{\link{psst}}.
  
  Model fitting: \code{\link{ppm}}.
}
\examples{
    data(cells)
    fit0 <- ppm(cells, ~1) # uniform Poisson
    G0 <- Gcom(fit0)
    G0
    plot(G0)
# uniform Poisson is clearly not correct

# Hanisch estimates only
    plot(Gcom(fit0), cbind(han, hcom) ~ r)

    fit1 <- ppm(cells, ~1, Strauss(0.08))
    plot(Gcom(fit1), cbind(han, hcom) ~ r)

# Try adjusting interaction distance

    fit2 <- update(fit1, Strauss(0.10))
    plot(Gcom(fit2), cbind(han, hcom) ~ r)

    G3 <- Gcom(cells, interaction=Strauss(0.12))
    plot(G3, cbind(han, hcom) ~ r)
}
\keyword{spatial}
\keyword{models}
