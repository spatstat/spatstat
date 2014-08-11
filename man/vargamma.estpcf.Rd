\name{vargamma.estpcf}
\alias{vargamma.estpcf}
\title{Fit the Neyman-Scott Cluster Point Process with Variance Gamma kernel}
\description{
  Fits the Neyman-Scott cluster point process, with Variance Gamma
  kernel, to a point pattern dataset by the Method of Minimum Contrast,
  using the pair correlation function.
}
\usage{
vargamma.estpcf(X, startpar=c(kappa=1,eta=1), nu.ker = -1/4, lambda=NULL,
            q = 1/4, p = 2, rmin = NULL, rmax = NULL, ..., pcfargs = list())
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
  \item{nu.ker}{
    Numerical value controlling the shape of the tail of the clusters.
    A number greater than \code{-1/2}.
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
  This algorithm fits the Neyman-Scott Cluster point process model
  with Variance Gamma kernel (Jalilian et al, 2011)
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

  The algorithm fits the Neyman-Scott Cluster point process
  with Variance Gamma kernel to \code{X},
  by finding the parameters of the model
  which give the closest match between the
  theoretical pair correlation function of the model
  and the observed pair correlation function.
  For a more detailed explanation of the Method of Minimum Contrast,
  see \code{\link{mincontrast}}.
  
  The Neyman-Scott cluster point process with Variance Gamma
  kernel is described in Jalilian et al (2011).
  It is a cluster process formed by taking a 
  pattern of parent points, generated according to a Poisson process
  with intensity \eqn{\kappa}{kappa}, and around each parent point,
  generating a random number of offspring points, such that the
  number of offspring of each parent is a Poisson random variable with mean
  \eqn{\mu}{mu}, and the locations of the offspring points of one parent
  have a common distribution described in Jalilian et al (2011).

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

  The corresponding model can be simulated using \code{\link{rVarGamma}}.
  
  The parameter \code{eta} appearing in \code{startpar} is equivalent to the
  scale parameter \code{omega} used in  \code{\link{rVarGamma}}.
  
  Homogeneous or inhomogeneous Neyman-Scott/VarGamma models can also be
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
  Jalilian, A., Guan, Y. and Waagepetersen, R. (2011)
  Decomposition of variance for spatial Cox processes.
  \emph{Scandinavian Journal of Statistics} \bold{40}, 119-137.

  Waagepetersen, R. (2007)
  An estimating function approach to inference for
  inhomogeneous Neyman-Scott processes.
  \emph{Biometrics} \bold{63}, 252--258.
}
\author{Abdollah Jalilian and Rasmus Waagepetersen.
  Adapted for \pkg{spatstat} by Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
}
\seealso{
  \code{\link{kppm}},
  \code{\link{vargamma.estK}},
  \code{\link{lgcp.estpcf}},
  \code{\link{thomas.estpcf}},
  \code{\link{cauchy.estpcf}},
  \code{\link{mincontrast}},
  \code{\link{pcf}},
  \code{\link{pcfmodel}}.

  \code{\link{rVarGamma}} to simulate the model.
}
\examples{
    u <- vargamma.estpcf(redwood)
    u
    plot(u, legendpos="topright")
}
\keyword{spatial}
\keyword{models}
