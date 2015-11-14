\name{cauchy.estK}
\alias{cauchy.estK}
\title{Fit the Neyman-Scott cluster process with Cauchy kernel}
\description{
  Fits the Neyman-Scott Cluster point process with Cauchy kernel
  to a point pattern dataset by the Method of Minimum Contrast.
}
\usage{
cauchy.estK(X, startpar=c(kappa=1,scale=1), lambda=NULL,
            q = 1/4, p = 2, rmin = NULL, rmax = NULL, ...)
}
\arguments{
  \item{X}{
    Data to which the model will be fitted.
    Either a point pattern or a summary statistic.
    See Details.
  }
  \item{startpar}{
    Vector of starting values for the parameters of the model.
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
}
\details{
  This algorithm fits the Neyman-Scott cluster point process model
  with Cauchy kernel to a point pattern dataset
  by the Method of Minimum Contrast, using the \eqn{K} function.

  The argument \code{X} can be either
  \describe{
    \item{a point pattern:}{An object of class \code{"ppp"}
      representing a point pattern dataset. 
      The \eqn{K} function of the point pattern will be computed
      using \code{\link{Kest}}, and the method of minimum contrast
      will be applied to this.
    }
    \item{a summary statistic:}{An object of class \code{"fv"} containing
      the values of a summary statistic, computed for a point pattern
      dataset. The summary statistic should be the \eqn{K} function,
      and this object should have been obtained by a call to
      \code{\link{Kest}} or one of its relatives.
    }
  }

  The algorithm fits the Neyman-Scott cluster point process
  with Cauchy kernel to \code{X},
  by finding the parameters of the Matern Cluster model
  which give the closest match between the
  theoretical \eqn{K} function of the Matern Cluster process
  and the observed \eqn{K} function.
  For a more detailed explanation of the Method of Minimum Contrast,
  see \code{\link{mincontrast}}.
  
  The model is described in Jalilian et al (2013).
  It is a cluster process formed by taking a 
  pattern of parent points, generated according to a Poisson process
  with intensity \eqn{\kappa}{kappa}, and around each parent point,
  generating a random number of offspring points, such that the
  number of offspring of each parent is a Poisson random variable with mean
  \eqn{\mu}{mu}, and the locations of the offspring points of one parent
  follow a common distribution described in Jalilian et al (2013).

  If the argument \code{lambda} is provided, then this is used
  as the value of the point process intensity \eqn{\lambda}{lambda}.
  Otherwise, if \code{X} is a
  point pattern, then  \eqn{\lambda}{lambda}
  will be estimated from \code{X}. 
  If \code{X} is a summary statistic and \code{lambda} is missing,
  then the intensity \eqn{\lambda}{lambda} cannot be estimated, and
  the parameter \eqn{\mu}{mu} will be returned as \code{NA}.

  The remaining arguments \code{rmin,rmax,q,p} control the
  method of minimum contrast; see \code{\link{mincontrast}}.

  The corresponding model can be simulated using \code{\link{rCauchy}}.

  For computational reasons, the optimisation procedure uses the parameter 
  \code{eta2}, which is equivalent to \code{4 * scale^2}
  where \code{scale} is the scale parameter for the model
  as used in \code{\link{rCauchy}}.
  
  Homogeneous or inhomogeneous Neyman-Scott/Cauchy models can also be
  fitted using the function \code{\link{kppm}} and the fitted models
  can be simulated using \code{\link{simulate.kppm}}.

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
  Ghorbani, M. (2012) Cauchy cluster process.
  \emph{Metrika}, to appear.

  Jalilian, A., Guan, Y. and Waagepetersen, R. (2013)
  Decomposition of variance for spatial Cox processes.
  \emph{Scandinavian Journal of Statistics} \bold{40}, 119-137.

  Waagepetersen, R. (2007)
  An estimating function approach to inference for
  inhomogeneous Neyman-Scott processes.
  \emph{Biometrics} \bold{63}, 252--258.
}
\author{Abdollah Jalilian and Rasmus Waagepetersen.
  Adapted for \pkg{spatstat} by \adrian
  
  
}
\seealso{
  \code{\link{kppm}},
  \code{\link{cauchy.estpcf}},
  \code{\link{lgcp.estK}},
  \code{\link{thomas.estK}},
  \code{\link{vargamma.estK}},
  \code{\link{mincontrast}},
  \code{\link{Kest}},
  \code{\link{Kmodel}}.

  \code{\link{rCauchy}} to simulate the model.
}
\examples{
    u <- cauchy.estK(redwood)
    u
    plot(u)
}
\keyword{spatial}
\keyword{models}
