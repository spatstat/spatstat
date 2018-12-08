\name{kppm}
\alias{kppm}
\alias{kppm.formula}
\alias{kppm.ppp}
\alias{kppm.quad}
\concept{point process model}
\concept{Cox point process}
\concept{cluster process}
\concept{Neyman-Scott cluster process}
\title{Fit Cluster or Cox Point Process Model}
\description{
  Fit a homogeneous or inhomogeneous cluster process or
  Cox point process model to a point pattern.
}
\usage{
  kppm(X, \dots)

  \method{kppm}{formula}(X,
                clusters = c("Thomas","MatClust","Cauchy","VarGamma","LGCP"),
                \dots,
                data=NULL)

  \method{kppm}{ppp}(X,
       trend = ~1,
       clusters = c("Thomas","MatClust","Cauchy","VarGamma","LGCP"),
       data = NULL,
       ...,
       covariates=data,
       subset,
       method = c("mincon", "clik2", "palm"),
       improve.type = c("none", "clik1", "wclik1", "quasi"),
       improve.args = list(),
       weightfun=NULL,
       control=list(),
       algorithm="Nelder-Mead",
       statistic="K",
       statargs=list(),
       rmax = NULL,
       covfunargs=NULL,
       use.gam=FALSE,
       nd=NULL, eps=NULL)

\method{kppm}{quad}(X,
       trend = ~1,
       clusters = c("Thomas","MatClust","Cauchy","VarGamma","LGCP"),
       data = NULL,
       ...,
       covariates=data,
       subset,
       method = c("mincon", "clik2", "palm"),
       improve.type = c("none", "clik1", "wclik1", "quasi"),
       improve.args = list(),
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
  \item{X}{
    A point pattern dataset (object of class \code{"ppp"} or
    \code{"quad"}) to which the model should be fitted, or a
    \code{formula} in the \R language defining the model. See Details.
  }
  \item{trend}{
    An \R formula, with no left hand side,
    specifying the form of the log intensity.
  }
  \item{clusters}{
    Character string determining the cluster model.
    Partially matched.
    Options are \code{"Thomas"}, \code{"MatClust"},
    \code{"Cauchy"}, \code{"VarGamma"} and \code{"LGCP"}.
  }
  \item{data,covariates}{
    The values of spatial covariates (other than the Cartesian
    coordinates) required by the model.
    A named list of pixel images, functions, windows,
    tessellations or numeric constants.
  }
  \item{\dots}{
    Additional arguments. See Details.
  }
  \item{subset}{
    Optional.
    A subset of the spatial domain,
    to which the model-fitting should be restricted.
    A window (object of class \code{"owin"})
    or a logical-valued pixel image (object of class \code{"im"}),
    or an expression (possibly involving the names of entries in \code{data})
    which can be evaluated to yield a window or pixel image.
  }
  \item{method}{
    The fitting method. Either 
    \code{"mincon"} for minimum contrast,
    \code{"clik2"} for second order composite likelihood,
    or \code{"palm"} for Palm likelihood.
    Partially matched.
  }
  \item{improve.type}{
    Method for updating the initial estimate of the trend.
    Initially the trend is estimated as if the process
    is an inhomogeneous Poisson process.
    The default, \code{improve.type = "none"}, is to use this initial estimate.
    Otherwise, the trend estimate is
    updated by \code{\link{improve.kppm}}, using information
    about the pair correlation function.
    Options are \code{"clik1"}
    (first order composite likelihood, essentially equivalent to \code{"none"}),
    \code{"wclik1"} (weighted first order composite likelihood) and
    \code{"quasi"} (quasi likelihood).
  }
  \item{improve.args}{
    Additional arguments passed to \code{\link{improve.kppm}} when
    \code{improve.type != "none"}. See Details.
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
  This function fits a clustered point process model to the
  point pattern dataset \code{X}. 

  The model may be either a \emph{Neyman-Scott cluster process}
  or another \emph{Cox process}.
  The type of model is determined by the argument \code{clusters}.
  Currently the options 
  are \code{clusters="Thomas"} for the Thomas process,
  \code{clusters="MatClust"} for the Matern cluster process,
  \code{clusters="Cauchy"} for the Neyman-Scott cluster process
  with Cauchy kernel,
  \code{clusters="VarGamma"} for the Neyman-Scott cluster process
  with Variance Gamma kernel (requires an additional argument \code{nu}
  to be passed through the dots; see \code{\link{rVarGamma}} for details),
  and \code{clusters="LGCP"} for the log-Gaussian Cox process (may
  require additional arguments passed through \code{\dots}; see
  \code{\link{rLGCP}} for details on argument names).
  The first four models are Neyman-Scott cluster processes.
  
  The algorithm first estimates the intensity function
  of the point process using \code{\link{ppm}}.
  The argument \code{X} may be a point pattern
  (object of class \code{"ppp"}) or a quadrature scheme
  (object of class \code{"quad"}). The intensity is specified by
  the \code{trend} argument.
  If the trend formula is \code{~1} (the default)
  then the model is \emph{homogeneous}. The algorithm begins by
  estimating the intensity as the number of points divided by
  the area of the window.
  Otherwise, the model is \emph{inhomogeneous}. 
  The algorithm begins by fitting a Poisson process with log intensity
  of the form specified by the formula \code{trend}.
  (See \code{\link{ppm}} for further explanation).

  The argument \code{X} may also be a \code{formula} in the
  \R language. The right hand side of the formula gives the
  \code{trend} as described above. The left hand side of the formula
  gives the point pattern dataset to which the model should be fitted.

  If \code{improve.type="none"} this is the final estimate of the
  intensity. Otherwise, the intensity estimate is updated, as explained in
  \code{\link{improve.kppm}}. Additional arguments to
  \code{\link{improve.kppm}} are passed as a named list in
  \code{improve.args}.
  
  The clustering parameters of the model are then fitted
  either by minimum contrast estimation, or by maximum
  composite likelihood.

  \describe{
   \item{Minimum contrast:}{
      If \code{method = "mincon"} (the default) clustering parameters of
      the model will be fitted
      by minimum contrast estimation, that is, by matching the theoretical
      \eqn{K}-function of the model to the empirical \eqn{K}-function
      of the data, as explained in \code{\link{mincontrast}}.

      For a homogeneous model (\code{ trend = ~1 })
      the empirical \eqn{K}-function of the data is computed
      using \code{\link{Kest}},
      and the parameters of the cluster model are estimated by
      the method of minimum contrast.

      For an inhomogeneous model, 
      the inhomogeneous \eqn{K} function is estimated
      by \code{\link{Kinhom}} using the fitted intensity.
      Then the parameters of the cluster model
      are estimated by the method of minimum contrast using the
      inhomogeneous \eqn{K} function. This two-step estimation
      procedure is due to Waagepetersen (2007).
  
      If \code{statistic="pcf"} then instead of using the
      \eqn{K}-function, the algorithm will use
      the pair correlation function \code{\link{pcf}} for homogeneous
      models and the inhomogeneous pair correlation function
      \code{\link{pcfinhom}} for inhomogeneous models.
      In this case, the smoothing parameters of the pair correlation
      can be controlled using the argument \code{statargs},
      as shown in the Examples.

      Additional arguments \code{\dots} will be passed to
      \code{\link{clusterfit}} to control the minimum contrast fitting
      algorithm.
    }
    \item{Composite likelihood:}{
      If \code{method = "clik2"} the clustering parameters of the
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
      If \code{method = "palm"} the clustering parameters of the
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

  Fitting the LGCP model requires the \pkg{RandomFields} package,
  except in the default case where the exponential covariance
  is assumed.
}
\section{Log-Gaussian Cox Models}{
  To fit a log-Gaussian Cox model with non-exponential covariance,
  specify \code{clusters="LGCP"} and use additional arguments
  to specify the covariance structure. These additional arguments can
  be given individually in the call to \code{kppm}, or they can be
  collected together in a list called \code{covmodel}.

  For example a Matern model with parameter \eqn{\nu=0.5} could be specified
  either by \code{kppm(X, clusters="LGCP", model="matern", nu=0.5)} or by
  \code{kppm(X, clusters="LGCP", covmodel=list(model="matern", nu=0.5))}.

  The argument \code{model} specifies the type of covariance
  model: the default is \code{model="exp"} for an exponential covariance.
  Alternatives include \code{"matern"}, \code{"cauchy"} and \code{"spheric"}.
  Model names correspond to functions beginning with \code{RM} in the
  \pkg{RandomFields} package: for example \code{model="matern"}
  corresponds to the function \code{RMmatern} in the 
  \pkg{RandomFields} package.
  
  Additional arguments are passed to the
  relevant function in the \pkg{RandomFields} package:
  for example if \code{model="matern"} then the additional argument
  \code{nu} is required, and is passed to the function
  \code{RMmatern} in the \pkg{RandomFields} package.

  Note that it is not possible to use \emph{anisotropic} covariance models
  because the \code{kppm} technique assumes the pair correlation function
  is isotropic.
}
\value{
  An object of class \code{"kppm"} representing the fitted model.
  There are methods for printing, plotting, predicting, simulating
  and updating objects of this class.
}
\section{Error and warning messages}{
  See \code{\link{ppm.ppp}} for a list of common error messages
  and warnings originating from the first stage of model-fitting.
}
\seealso{
  Methods for \code{kppm} objects:
  \code{\link{plot.kppm}},
  \code{\link{fitted.kppm}},
  \code{\link{predict.kppm}},
  \code{\link{simulate.kppm}},
  \code{\link{update.kppm}},
  \code{\link{vcov.kppm}},
  \code{\link[spatstat:methods.kppm]{methods.kppm}},
  \code{\link{as.ppm.kppm}},
  \code{\link{as.fv.kppm}},
  \code{\link{Kmodel.kppm}},
  \code{\link{pcfmodel.kppm}}.

  Minimum contrast fitting algorithm:
  higher level interface \code{\link{clusterfit}};
  low-level algorithm \code{\link{mincontrast}}.

  Alternative fitting algorithms:
  \code{\link{thomas.estK}},
  \code{\link{matclust.estK}},
  \code{\link{lgcp.estK}},
  \code{\link{cauchy.estK}},
  \code{\link{vargamma.estK}},
  \code{\link{thomas.estpcf}},
  \code{\link{matclust.estpcf}},
  \code{\link{lgcp.estpcf}},
  \code{\link{cauchy.estpcf}},
  \code{\link{vargamma.estpcf}},

  Summary statistics:
  \code{\link{Kest}},
  \code{\link{Kinhom}},
  \code{\link{pcf}},
  \code{\link{pcfinhom}}.

  See also \code{\link{ppm}}
}
\references{
  Guan, Y. (2006) 
  A composite likelihood approach in fitting spatial point process
  models.
  \emph{Journal of the American Statistical Association}
  \bold{101}, 1502--1512.

  Jalilian, A., Guan, Y. and Waagepetersen, R. (2012)
  Decomposition of variance for spatial Cox processes.
  \emph{Scandinavian Journal of Statistics} \bold{40}, 119--137.

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
  # method for point patterns
  kppm(redwood, ~1, "Thomas")
  # method for formulas
  kppm(redwood ~ 1, "Thomas")

  kppm(redwood ~ 1, "Thomas", method="c")
  kppm(redwood ~ 1, "Thomas", method="p")

  kppm(redwood ~ x, "MatClust") 
  kppm(redwood ~ x, "MatClust", statistic="pcf", statargs=list(stoyan=0.2)) 
  kppm(redwood ~ x, cluster="Cauchy", statistic="K")
  kppm(redwood, cluster="VarGamma", nu = 0.5, statistic="pcf")

  # LGCP models
  kppm(redwood ~ 1, "LGCP", statistic="pcf")
  if(require("RandomFields")) {
    kppm(redwood ~ x, "LGCP", statistic="pcf",
                              model="matern", nu=0.3,
                              control=list(maxit=10))
  }

  # fit with composite likelihood method
  kppm(redwood ~ x, "VarGamma", method="clik2", nu.ker=-3/8)

  # fit intensity with quasi-likelihood method
  kppm(redwood ~ x, "Thomas", improve.type = "quasi")
}
\author{
  \spatstatAuthors,
  with contributions from Abdollah Jalilian and Rasmus Waagepetersen.
}
\keyword{spatial}
\keyword{models}

