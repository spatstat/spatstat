\name{lgcp.estpcf}
\alias{lgcp.estpcf}
\title{Fit a Log-Gaussian Cox Point Process by Minimum Contrast}
\description{
  Fits a log-Gaussian Cox point process model
  (with exponential covariance function)
  to a point pattern dataset by the Method of Minimum Contrast
  using the pair correlation function.
}
\usage{
lgcp.estpcf(X,
            startpar=c(sigma2=1,alpha=1),
            covmodel=list(model="exponential"),
            lambda=NULL,
            q = 1/4, p = 2, rmin = NULL, rmax = NULL, ..., pcfargs=list())
}
\arguments{
  \item{X}{
    Data to which the model will be fitted.
    Either a point pattern or a summary statistic.
    See Details.
  }
  \item{startpar}{
    Vector of starting values for the parameters of the
    log-Gaussian Cox process model.
  }
  \item{covmodel}{
    Specification of the covariance model
    for the log-Gaussian field. See Details.
  }
  \item{lambda}{
    Optional. An estimate of the intensity of the point process.
  }
  \item{q,p}{
    Optional. Exponents for the contrast criterion.
  }
  \item{rmin, rmax}{
    Optional. The interval of \eqn{r} values for the contrast criterion.
  }
  \item{\dots}{
    Optional arguments passed to \code{\link[stats]{optim}}
    to control the optimisation algorithm. See Details.
  }
  \item{pcfargs}{
    Optional list containing arguments passed to \code{\link{pcf.ppp}}
    to control the smoothing in the estimation of the
    pair correlation function.
  }
}
\details{
  This algorithm fits a log-Gaussian Cox point process model
  to a point pattern dataset
  by the Method of Minimum Contrast, using the pair correlation function.

  The argument \code{X} can be either
  \describe{
    \item{a point pattern:}{An object of class \code{"ppp"}
      representing a point pattern dataset. 
      The pair correlation function of the point pattern will be computed
      using \code{\link{pcf}}, and the method of minimum contrast
      will be applied to this.
    }
    \item{a summary statistic:}{An object of class \code{"fv"} containing
      the values of a summary statistic, computed for a point pattern
      dataset. The summary statistic should be the pair correlation function,
      and this object should have been obtained by a call to
      \code{\link{pcf}} or one of its relatives.
    }
  }
  The algorithm fits a log-Gaussian Cox point process (LGCP)
  model to \code{X},  by finding the parameters of the LGCP model
  which give the closest match between the
  theoretical pair correlation function of the LGCP model
  and the observed pair correlation function.
  For a more detailed explanation of the Method of Minimum Contrast,
  see \code{\link{mincontrast}}.

  The model fitted is a stationary, isotropic log-Gaussian Cox process
  (\ifelse{latex}{\out{M\o ller}}{Moller} and Waagepetersen, 2003, pp. 72-76).
  To define this process we start with
  a stationary Gaussian random field \eqn{Z} in the two-dimensional plane,
  with constant mean \eqn{\mu}{mu} and covariance function \eqn{C(r)}.
  Given \eqn{Z}, we generate a Poisson point process \eqn{Y} with intensity
  function \eqn{\lambda(u) = \exp(Z(u))}{lambda(u) = exp(Z(u))} at
  location \eqn{u}. Then \eqn{Y} is a log-Gaussian Cox process.

  The theoretical pair correlation function of the LGCP is
  \deqn{
    g(r) = \exp(C(s))
  }{
    g(r) = exp(C(s))
  }
  The intensity of the LGCP is 
  \deqn{
    \lambda = \exp(\mu + \frac{C(0)}{2}).
  }{
    lambda= exp(mu + C(0)/2).
  }
  
  The covariance function \eqn{C(r)} takes the form
  \deqn{
    C(r) = \sigma^2 c(r/\alpha)
  }{
    C(r) = sigma^2 * c(-r/alpha)
  }
  where \eqn{\sigma^2}{sigma^2} and \eqn{\alpha}{alpha} are parameters
  controlling the strength and the scale of autocorrelation,
  respectively, and \eqn{c(r)} is a known covariance function
  determining the shape of the covariance. 
  The strength and scale parameters
  \eqn{\sigma^2}{sigma^2} and \eqn{\alpha}{alpha}
  will be estimated by the algorithm.
  The template covariance function \eqn{c(r)} must be specified
  as explained below.
  
  In this algorithm, the Method of Minimum Contrast is first used to find
  optimal values of the parameters \eqn{\sigma^2}{sigma^2}
  and \eqn{\alpha}{alpha^2}. Then the remaining parameter
  \eqn{\mu}{mu} is inferred from the estimated intensity
  \eqn{\lambda}{lambda}.

  The template covariance function \eqn{c(r)} is specified
  using the argument \code{covmodel}. This should be of the form
  \code{list(model="modelname", \dots)} where
  \code{modelname} is a string identifying the template model
  as explained below, and  \code{\dots} are optional arguments of the
  form \code{tag=value} giving the values of parameters controlling the
  \emph{shape} of the template model.
  The default is the exponential covariance
  \eqn{c(r) = e^{-r}}{c(r) = e^(-r)}
  so that the scaled covariance is 
  \deqn{
    C(r) = \sigma^2 e^{-r/\alpha}.
  }{
    C(r) = sigma^2 * exp(-r/alpha).
  }
  To determine the template model, the string \code{"modelname"} will be
  prefixed by \code{"RM"} and the code will search for
  a function of this name in the \pkg{RandomFields} package.
  For a list of available models see 
  \code{\link[RandomFields]{RMmodel}} in the
  \pkg{RandomFields} package. For example the
  Matern covariance with exponent \eqn{\nu=0.3}{nu = 0.3} is specified
  by \code{covmodel=list(model="matern", nu=0.3)} corresponding
  to the function \code{RMmatern} in the \pkg{RandomFields} package.
  
  If the argument \code{lambda} is provided, then this is used
  as the value of \eqn{\lambda}{lambda}. Otherwise, if \code{X} is a
  point pattern, then  \eqn{\lambda}{lambda}
  will be estimated from \code{X}. 
  If \code{X} is a summary statistic and \code{lambda} is missing,
  then the intensity \eqn{\lambda}{lambda} cannot be estimated, and
  the parameter \eqn{\mu}{mu} will be returned as \code{NA}.

  The remaining arguments \code{rmin,rmax,q,p} control the
  method of minimum contrast; see \code{\link{mincontrast}}.

  The optimisation algorithm can be controlled through the
  additional arguments \code{"..."} which are passed to the
  optimisation function \code{\link[stats]{optim}}. For example,
  to constrain the parameter values to a certain range,
  use the argument \code{method="L-BFGS-B"} to select an optimisation
  algorithm that respects box constraints, and use the arguments
  \code{lower} and \code{upper} to specify (vectors of) minimum and
  maximum values for each parameter.
}
\value{
  An object of class \code{"minconfit"}. There are methods for printing
  and plotting this object. It contains the following main components:
  \item{par }{Vector of fitted parameter values.}
  \item{fit }{Function value table (object of class \code{"fv"})
    containing the observed values of the summary statistic
    (\code{observed}) and the theoretical values of the summary
    statistic computed from the fitted model parameters.
  }
}
\references{
  \ifelse{latex}{\out{M\o ller}}{Moller}, J, Syversveen, A. and Waagepetersen, R. (1998)
  Log Gaussian Cox Processes.
  \emph{Scandinavian Journal of Statistics} \bold{25}, 451--482.
  
  \ifelse{latex}{\out{M\o ller}}{Moller}, J. and Waagepetersen, R. (2003).
  Statistical Inference and Simulation for Spatial Point Processes.
  Chapman and Hall/CRC, Boca Raton.

  Waagepetersen, R. (2007)
  An estimating function approach to inference for
  inhomogeneous Neyman-Scott processes.
  \emph{Biometrics} \bold{63}, 252--258.
}
\author{
  Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/},
  with modifications by Shen Guochun and 
  Rasmus Waagepetersen
  \email{rw@math.auc.dk}
}
\seealso{
  \code{\link{lgcp.estK}} for alternative method of fitting LGCP.
  
  \code{\link{matclust.estpcf}},
  \code{\link{thomas.estpcf}} for other models.
  
  \code{\link{mincontrast}} for the generic minimum contrast
  fitting algorithm, including important parameters that affect
  the accuracy of the fit.
  
  \code{\link[RandomFields]{RMmodel}} in the
  \pkg{RandomFields} package, for covariance function models.
  
  \code{\link{pcf}} for the pair correlation function.
}
\examples{
    data(redwood)
    u <- lgcp.estpcf(redwood, c(sigma2=0.1, alpha=1))
    u
    plot(u)
    if(require(RandomFields)) {
      lgcp.estpcf(redwood, covmodel=list(model="matern", nu=0.3))
    }
}
\keyword{spatial}
\keyword{models}
