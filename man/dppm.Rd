\name{dppm}
\alias{dppm}
\title{Fit Determinantal Point Process Model}
\description{
  Fit a determinantal point process model to a point pattern.
}
\usage{
  dppm(formula, family, data=NULL,
       ...,
       startpar = NULL,
       method = c("mincon", "clik2", "palm"),
       weightfun=NULL,
       control=list(),
       algorithm="Nelder-Mead",
       statistic="K",
       statargs=list(),
       rmax = NULL,
       covfunargs=NULL,
       use.gam=FALSE,
       nd=NULL, eps=NULL)
}
\arguments{
  \item{formula}{
    A \code{formula} in the \R language
    specifying the data (on the left side) and the
    form of the model to be fitted (on the right side).
    For a stationary model it suffices to provide a point pattern
    without a formula. See Details.
  }
  \item{family}{
    Information specifying the family of point processes
    to be used in the model.
    Typically one of the family functions
    \code{\link{dppGauss}}, \code{\link{dppMatern}},
    \code{\link{dppCauchy}}, \code{\link{dppBessel}}
    or \code{\link{dppPowerExp}}.
    Alternatively a character string giving the name
    of a family function, or the result of calling one of the
    family functions. See Details.
  }
  \item{data}{
    The values of spatial covariates (other than the Cartesian
    coordinates) required by the model.
    A named list of pixel images, functions, windows,
    tessellations or numeric constants.
  }
  \item{\dots}{
    Additional arguments. See Details.
  }
  \item{startpar}{
    Named vector of starting parameter values for the optimization.
  }
  \item{method}{
    The fitting method. Either
    \code{"mincon"} for minimum contrast,
    \code{"clik2"} for second order composite likelihood,
    or \code{"palm"} for Palm likelihood.
    Partially matched.
  }
  \item{weightfun}{
    Optional weighting function \eqn{w}
    in the composite likelihood or Palm likelihood.
    A \code{function} in the \R language.
    See Details.
  }
  \item{control}{
    List of control parameters passed to the optimization function
    \code{\link[stats]{optim}}.
  }
  \item{algorithm}{
    Character string determining the mathematical optimisation algorithm
    to be used by \code{\link[stats]{optim}}. See
    the argument \code{method} of \code{\link[stats]{optim}}.
  }
  \item{statistic}{
    Name of the summary statistic to be used
    for minimum contrast estimation: either \code{"K"} or \code{"pcf"}.
  }
  \item{statargs}{
    Optional list of arguments to be used when calculating
    the \code{statistic}. See Details.
  }
  \item{rmax}{
    Maximum value of interpoint distance
    to use in the composite likelihood.
  }
  \item{covfunargs,use.gam,nd,eps}{
    Arguments passed to \code{\link{ppm}} when fitting the intensity.
  }
}
\details{
  This function fits a determinantal point process model to a
  point pattern dataset as described in Lavancier et al. (2015).

  The model to be fitted is specified by the arguments
  \code{formula} and \code{family}.
  
  The argument \code{formula} should normally be a \code{formula} in the
  \R language. The left hand side of the formula
  specifies the point pattern dataset to which the model should be fitted.
  This should be a single argument which may be a point pattern
  (object of class \code{"ppp"}) or a quadrature scheme
  (object of class \code{"quad"}). The right hand side of the formula is called
  the \code{trend} and specifies the form of the
  \emph{logarithm of the intensity} of the process.
  Alternatively the argument \code{formula} may be a point pattern or quadrature
  scheme, and the trend formula is taken to be \code{~1}.

  The argument \code{family} specifies the family of point processes
  to be used in the model.
  It is typically one of the family functions
  \code{\link{dppGauss}}, \code{\link{dppMatern}},
  \code{\link{dppCauchy}}, \code{\link{dppBessel}}
  or \code{\link{dppPowerExp}}. 
  Alternatively it may be a character string giving the name
  of a family function, or the result of calling one of the
  family functions. A family function belongs to class
  \code{"detpointprocfamilyfun"}. The result of calling a family
  function is a point process family, which belongs to class
  \code{"detpointprocfamily"}.   
  
  The algorithm first estimates the intensity function
  of the point process using \code{\link{ppm}}.
  If the trend formula is \code{~1}
  (the default if a point pattern or quadrature
  scheme is given rather than a \code{"formula"})
  then the model is \emph{homogeneous}. The algorithm begins by
  estimating the intensity as the number of points divided by
  the area of the window.
  Otherwise, the model is \emph{inhomogeneous}.
  The algorithm begins by fitting a Poisson process with log intensity
  of the form specified by the formula \code{trend}.
  (See \code{\link{ppm}} for further explanation).

  The interaction parameters of the model are then fitted
  either by minimum contrast estimation, or by maximum
  composite likelihood.

  \describe{
   \item{Minimum contrast:}{
      If \code{method = "mincon"} (the default) interaction parameters of
      the model will be fitted
      by minimum contrast estimation, that is, by matching the theoretical
      \eqn{K}-function of the model to the empirical \eqn{K}-function
      of the data, as explained in \code{\link{mincontrast}}.

      For a homogeneous model (\code{ trend = ~1 })
      the empirical \eqn{K}-function of the data is computed
      using \code{\link{Kest}},
      and the interaction parameters of the model are estimated by
      the method of minimum contrast.

      For an inhomogeneous model,
      the inhomogeneous \eqn{K} function is estimated
      by \code{\link{Kinhom}} using the fitted intensity.
      Then the interaction parameters of the model
      are estimated by the method of minimum contrast using the
      inhomogeneous \eqn{K} function. This two-step estimation
      procedure is heavily inspired by Waagepetersen (2007).

      If \code{statistic="pcf"} then instead of using the
      \eqn{K}-function, the algorithm will use
      the pair correlation function \code{\link{pcf}} for homogeneous
      models and the inhomogeneous pair correlation function
      \code{\link{pcfinhom}} for inhomogeneous models.
      In this case, the smoothing parameters of the pair correlation
      can be controlled using the argument \code{statargs},
      as shown in the Examples.

      Additional arguments \code{\dots} will be passed to
      \code{\link{mincontrast}} to control the minimum contrast fitting
      algorithm.
    }
    \item{Composite likelihood:}{
      If \code{method = "clik2"} the interaction parameters of the
      model will be fitted by maximising the second-order composite likelihood
      (Guan, 2006). The log composite likelihood is
      \deqn{
	\sum_{i,j} w(d_{ij}) \log\rho(d_{ij}; \theta)
	- \left( \sum_{i,j} w(d_{ij}) \right)
	\log \int_D \int_D w(\|u-v\|) \rho(\|u-v\|; \theta)\, du\, dv
      }{
	sum[i,j] w(d[i,j]) log(rho(d[i,j]; theta))
	- (sum[i,j] w(d[i,j]))
	log(integral[D,D] w(||u-v||) rho(||u-v||; theta) du dv)
      }
      where the sums are taken over all pairs of data points
      \eqn{x_i, x_j}{x[i], x[j]} separated by a distance
      \eqn{d_{ij} = \| x_i - x_j\|}{d[i,j] = ||x[i] - x[j]||}
      less than \code{rmax},
      and the double integral is taken over all pairs of locations
      \eqn{u,v} in the spatial window of the data.
      Here \eqn{\rho(d;\theta)}{rho(d;theta)} is the
      pair correlation function of the model with
      cluster parameters \eqn{\theta}{theta}.

      The function \eqn{w} in the composite likelihood
      is a weighting function and may be chosen arbitrarily.
      It is specified by the argument \code{weightfun}.
      If this is missing or \code{NULL} then the default is
      a threshold weight function,
      \eqn{w(d) = 1(d \le R)}{w(d) = 1(d <= R)}, where \eqn{R} is \code{rmax/2}.
    }
    \item{Palm likelihood:}{
      If \code{method = "palm"} the interaction parameters of the
      model will be fitted by maximising the Palm loglikelihood
      (Tanaka et al, 2008)
      \deqn{
	\sum_{i,j} w(x_i, x_j) \log \lambda_P(x_j \mid x_i; \theta)
	- \int_D w(x_i, u) \lambda_P(u \mid x_i; \theta) {\rm d} u
      }{
	sum[i,j] w(x[i], x[j]) log(lambdaP(x[j] | x[i]; theta)
	- integral[D] w(x[i], u) lambdaP(u | x[i]; theta) du
      }
      with the same notation as above. Here
      \eqn{\lambda_P(u|v;\theta}{lambdaP(u|v;theta)} is the Palm intensity of
      the model at location \eqn{u} given there is a point at \eqn{v}.
    }
  }
  In all three methods, the optimisation is performed by the generic
  optimisation algorithm \code{\link[stats]{optim}}.
  The behaviour of this algorithm can be modified using the
  argument \code{control}.
  Useful control arguments include
  \code{trace}, \code{maxit} and \code{abstol}
  (documented in the help for \code{\link[stats]{optim}}).

  Finally, it is also possible to fix any parameters desired before the
  optimisation by specifying them as \code{name=value}
  in the call to the family function. See Examples.
}
\value{
  An object of class \code{"dppm"} representing the fitted model.
  There are methods for printing, plotting, predicting and simulating
  objects of this class.
}
\seealso{
  methods for \code{dppm} objects:
  \code{\link{plot.dppm}},
  \code{\link{fitted.dppm}},
  \code{\link{predict.dppm}},
  \code{\link{simulate.dppm}},
  \code{\link{methods.dppm}},
  \code{\link{as.ppm.dppm}},
  \code{\link{Kmodel.dppm}},
  \code{\link{pcfmodel.dppm}}.

  Minimum contrast fitting algorithm:
  \code{\link{mincontrast}}.

  Deterimantal point process models:
  \code{\link{dppGauss}},
  \code{\link{dppMatern}},
  \code{\link{dppCauchy}},
  \code{\link{dppBessel}},
  \code{\link{dppPowerExp}},

  Summary statistics:
  \code{\link{Kest}},
  \code{\link{Kinhom}},
  \code{\link{pcf}},
  \code{\link{pcfinhom}}.

  See also \code{\link{ppm}}
}
\references{
  Lavancier, F. \ifelse{latex}{\out{M\o ller}}{Moller}, J. and Rubak, E. (2015)
  Determinantal point process models and statistical inference
  \emph{Journal of the Royal Statistical Society, Series B}
  \bold{77}, 853--977.

  Guan, Y. (2006)
  A composite likelihood approach in fitting spatial point process
  models.
  \emph{Journal of the American Statistical Association}
  \bold{101}, 1502--1512.

  Tanaka, U. and Ogata, Y. and Stoyan, D. (2008)
  Parameter estimation and model selection for
  Neyman-Scott point processes.
  \emph{Biometrical Journal} \bold{50}, 43--57.

  Waagepetersen, R. (2007)
  An estimating function approach to inference for
  inhomogeneous Neyman-Scott processes.
  \emph{Biometrics} \bold{63}, 252--258.
}
\examples{
  jpines <- residualspaper$Fig1

  dppm(jpines ~ 1, dppGauss)

  dppm(jpines ~ 1, dppGauss, method="c")
  dppm(jpines ~ 1, dppGauss, method="p")

  # Fixing the intensity to lambda=2 rather than the Poisson MLE 2.04:
  dppm(jpines ~ 1, dppGauss(lambda=2))

  if(interactive()) {
   # The following is quite slow (using K-function)
   dppm(jpines ~ x, dppMatern)
  }

   # much faster using pair correlation function
  dppm(jpines ~ x, dppMatern, statistic="pcf", statargs=list(stoyan=0.2))

  # Fixing the Matern shape parameter to nu=2 rather than estimating it:
  dppm(jpines ~ x, dppMatern(nu=2))
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/},
  Rolf Turner
  \email{r.turner@auckland.ac.nz}
  and Ege Rubak
  \email{rubak@math.aau.dk}
}
\keyword{spatial}
\keyword{models}

