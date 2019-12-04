\name{matclust.estK}
\alias{matclust.estK}
\title{Fit the Matern Cluster Point Process by Minimum Contrast}
\description{
  Fits the \Matern Cluster point process to a point pattern dataset
  by the Method of Minimum Contrast.
}
\usage{
matclust.estK(X, startpar=c(kappa=1,scale=1), lambda=NULL,
            q = 1/4, p = 2, rmin = NULL, rmax = NULL, ...)
}
\arguments{
  \item{X}{
    Data to which the \Matern Cluster model will be fitted.
    Either a point pattern or a summary statistic.
    See Details.
  }
  \item{startpar}{
    Vector of starting values for the parameters of the
    \Matern Cluster process.
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
  This algorithm fits the \Matern Cluster point process model
  to a point pattern dataset
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

  The algorithm fits the \Matern Cluster point process to \code{X},
  by finding the parameters of the \Matern Cluster model
  which give the closest match between the
  theoretical \eqn{K} function of the \Matern Cluster process
  and the observed \eqn{K} function.
  For a more detailed explanation of the Method of Minimum Contrast,
  see \code{\link{mincontrast}}.
  
  The \Matern Cluster point process is described in \Moller and Waagepetersen
  (2003, p. 62). It is a cluster process formed by taking a 
  pattern of parent points, generated according to a Poisson process
  with intensity \eqn{\kappa}{kappa}, and around each parent point,
  generating a random number of offspring points, such that the
  number of offspring of each parent is a Poisson random variable with mean
  \eqn{\mu}{mu}, and the locations of the offspring points of one parent
  are independent and uniformly distributed inside a circle of radius
  \eqn{R}{R} centred on the parent point, where \eqn{R}{R} is equal to
  the parameter \code{scale}. The named vector of stating values can use
  either \code{R} or \code{scale} as the name of the second component,
  but the latter is recommended for consistency with other cluster models.

  The theoretical \eqn{K}-function of the \Matern Cluster process is
  \deqn{
    K(r) = \pi r^2 + \frac 1 \kappa h(\frac{r}{2R})
  }{
    K(r) = pi r^2 + h(r/(2*R))/kappa
  }
  where the radius R is the parameter \code{scale} and
  \deqn{
    h(z) = 2 + \frac 1 \pi [ ( 8 z^2 - 4 ) \mbox{arccos}(z)
    - 2 \mbox{arcsin}(z)
    + 4 z \sqrt{(1 - z^2)^3}
    - 6 z \sqrt{1 - z^2}
    ]
  }{
    h(z) = 2 + (1/pi) * ((8 * z^2 - 4) * arccos(z) - 2 * arcsin(z) 
    + 4 * z * sqrt((1 - z^2)^3) - 6 * z * sqrt(1 - z^2))
  }
  for \eqn{z <= 1}, and \eqn{h(z) = 1} for \eqn{z > 1}.
  The theoretical intensity
  of the \Matern Cluster process
  is \eqn{\lambda = \kappa \mu}{lambda=kappa* mu}.

  In this algorithm, the Method of Minimum Contrast is first used to find
  optimal values of the parameters \eqn{\kappa}{kappa}
  and \eqn{R}{R}. Then the remaining parameter
  \eqn{\mu}{mu} is inferred from the estimated intensity
  \eqn{\lambda}{lambda}.

  If the argument \code{lambda} is provided, then this is used
  as the value of \eqn{\lambda}{lambda}. Otherwise, if \code{X} is a
  point pattern, then  \eqn{\lambda}{lambda}
  will be estimated from \code{X}. 
  If \code{X} is a summary statistic and \code{lambda} is missing,
  then the intensity \eqn{\lambda}{lambda} cannot be estimated, and
  the parameter \eqn{\mu}{mu} will be returned as \code{NA}.

  The remaining arguments \code{rmin,rmax,q,p} control the
  method of minimum contrast; see \code{\link{mincontrast}}.

  The \Matern Cluster process can be simulated, using
  \code{\link{rMatClust}}.

  Homogeneous or inhomogeneous \Matern Cluster models can also be
  fitted using the function \code{\link{kppm}}.

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
  \Moller, J. and Waagepetersen, R. (2003).
  Statistical Inference and Simulation for Spatial Point Processes.
  Chapman and Hall/CRC, Boca Raton.

  Waagepetersen, R. (2007)
  An estimating function approach to inference for
  inhomogeneous Neyman-Scott processes.
  \emph{Biometrics} \bold{63}, 252--258.
}
\author{Rasmus Waagepetersen
  \email{rw@math.auc.dk}
  Adapted for \pkg{spatstat} by \adrian
  
  
}
\seealso{
  \code{\link{kppm}},
  \code{\link{lgcp.estK}},
  \code{\link{thomas.estK}},
  \code{\link{mincontrast}},
  \code{\link{Kest}},
  \code{\link{rMatClust}} to simulate the fitted model.
}
\examples{
    data(redwood)
    u <- matclust.estK(redwood, c(kappa=10, scale=0.1))
    u
    plot(u)
}
\keyword{spatial}
\keyword{models}
