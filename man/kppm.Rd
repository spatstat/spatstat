\name{kppm}
\alias{kppm}
\title{Fit Cluster or Cox Point Process Model}
\description{
  Fit a homogeneous or inhomogeneous cluster process or
  Cox point process model to a point pattern.
}
\usage{
  kppm(X, trend = ~1,
       clusters = c("Thomas","MatClust","Cauchy","VarGamma","LGCP"),
       covariates = NULL,
       ...,
       method = c("mincon", "clik"),
       weightfun=NULL,
       control=list(),
       statistic="K",
       statargs=list(),
       rmax = NULL,
       covfunargs=NULL,
       use.gam=FALSE,
       nd=NULL, eps=NULL)
}
\arguments{
  \item{X}{
    Point pattern dataset (object of class \code{"ppp"}) to which the model
    should be fitted. See Details.
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
  \item{covariates}{
    The values of any spatial covariates (other than the Cartesian
    coordinates) required by the model.
    A named list of pixel images, functions, windows or numeric constants.
  }
  \item{\dots}{
    Additional arguments. See Details.
  }
  \item{method}{
    The fitting method. Either 
    \code{"mincon"} for minimum contrast,
    or \code{"clik"} for composite likelihood.
    Partially matched.
  }
  \item{weightfun}{
    Optional weighting function \eqn{w}
    in the composite likelihood. A \code{function} in the \R language.
    See Details.
  }
  \item{control}{
    List of control parameters passed to the optimization algorithm
    \code{\link[stats]{optim}}.
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

  The model may be either a \emph{Poisson cluster process}
  or another \emph{Cox process}.
  The type of model is determined by the argument \code{clusters}.
  Currently the options 
  are \code{clusters="Thomas"} for the Thomas process,
  \code{clusters="MatClust"} for the Matern cluster process,
  \code{clusters="Cauchy"} for the Neyman-Scott cluster process
  with Cauchy kernel,
  \code{clusters="VarGamma"} for the Neyman-Scott cluster process
  with Variance Gamma kernel,
  and \code{clusters="LGCP"} for the log-Gaussian Cox process.
  The first four models are Poisson cluster processes.

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
      the appropriate fitting function 
      \code{\link{thomas.estK}} or
      \code{\link{thomas.estpcf}} or
      \code{\link{matclust.estK}} or
      \code{\link{matclust.estpcf}} or
      \code{\link{lgcp.estK}} or
      \code{\link{lgcp.estpcf}} or
      \code{\link{cauchy.estK}} or
      \code{\link{cauchy.estpcf}} or
      \code{\link{vargamma.estK}} or
      \code{\link{vargamma.estpcf}}
      to control the minimum contrast fitting algorithm.
    }
    \item{Composite likelihood:}{
      If \code{method = "clik"} the clustering parameters of the
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
      If this is missing or \code{NULL} then \eqn{w} is taken to
      equal 1.
    }
  }
  In both methods, the optimisation is performed by the generic
  optimisation algorithm \code{\link[stats]{optim}}.
  The behaviour of this algorithm can be modified using the
  argument \code{control}.
  Useful control arguments include
  \code{trace}, \code{maxit} and \code{abstol}
  (documented in the help for \code{\link[stats]{optim}}).
}
\value{
  An object of class \code{"kppm"} representing the fitted model.
  There are methods for printing, plotting, predicting, simulating
  and updating objects of this class.
}
\seealso{
  methods for \code{kppm} objects:
  \code{\link{plot.kppm}},
  \code{\link{fitted.kppm}},
  \code{\link{predict.kppm}},
  \code{\link{simulate.kppm}},
  \code{\link{update.kppm}},
  \code{\link{vcov.kppm}},
  \code{\link[spatstat:methods.kppm]{methods.kppm}},
  \code{\link{as.ppm.kppm}},
  \code{\link{Kmodel.kppm}},
  \code{\link{pcfmodel.kppm}}.

  Minimum contrast fitting algorithms:
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
  \code{\link{mincontrast}}.

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
  \emph{Scandinavian Journal of Statistics}, in press.

  Waagepetersen, R. (2007)
  An estimating function approach to inference for
  inhomogeneous Neyman-Scott processes.
  \emph{Biometrics} \bold{63}, 252--258.
}
\examples{
  data(redwood)
  kppm(redwood, ~1, "Thomas")
  kppm(redwood, ~1, "Thomas", method="c")
  kppm(redwood, ~x, "MatClust") 
  kppm(redwood, ~x, "MatClust", statistic="pcf", statargs=list(stoyan=0.2)) 
  kppm(redwood, ~1, "LGCP", statistic="pcf")
  kppm(redwood, ~x, cluster="Cauchy", statistic="K")
  kppm(redwood, cluster="VarGamma", nu.ker = 0.5, statistic="pcf")
  if(require(RandomFields)) {
     kppm(redwood, ~x, "LGCP", statistic="pcf",
           covmodel=list(model="matern", nu=0.3),
           control=list(maxit=10))
  }
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/},
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
  with contributions from Abdollah Jalilian and Rasmus Waagepetersen.
}
\keyword{spatial}
\keyword{models}

