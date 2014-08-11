\name{thomas.estpcf}
\alias{thomas.estpcf}
\title{Fit the Thomas Point Process by Minimum Contrast}
\description{
  Fits the Thomas point process to a point pattern dataset by the Method of
  Minimum Contrast using the pair correlation function.
}
\usage{
thomas.estpcf(X, startpar=c(kappa=1,sigma2=1), lambda=NULL,
            q = 1/4, p = 2, rmin = NULL, rmax = NULL, ..., pcfargs=list())
}
\arguments{
  \item{X}{
    Data to which the Thomas model will be fitted.
    Either a point pattern or a summary statistic.
    See Details.
  }
  \item{startpar}{
    Vector of starting values for the parameters of the
    Thomas process.
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
  This algorithm fits the Thomas point process model to a point pattern dataset
  by the Method of Minimum Contrast, using the pair correlation function
  \code{\link{pcf}}.

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

  The algorithm fits the Thomas point process to \code{X},
  by finding the parameters of the Thomas model
  which give the closest match between the
  theoretical pair correlation function of the Thomas process
  and the observed pair correlation function.
  For a more detailed explanation of the Method of Minimum Contrast,
  see \code{\link{mincontrast}}.
  
  The Thomas point process is described in Moller and Waagepetersen
  (2003, pp. 61--62). It is a cluster process formed by taking a 
  pattern of parent points, generated according to a Poisson process
  with intensity \eqn{\kappa}{kappa}, and around each parent point,
  generating a random number of offspring points, such that the
  number of offspring of each parent is a Poisson random variable with mean
  \eqn{\mu}{mu}, and the locations of the offspring points of one parent
  are independent and isotropically Normally distributed around the parent
  point with standard deviation \eqn{\sigma}{sigma}.

  The theoretical pair correlation function of the Thomas process is
  \deqn{
    g(r) = 1 + \frac 1 {4\pi \kappa \sigma^2} \exp(-\frac{r^2}{4\sigma^2})).
  }{
    g(r) = 1 + exp(-r^2/(4 * sigma^2)))/(4 * pi * kappa * sigma^2).
  }
  The theoretical intensity
  of the Thomas process is \eqn{\lambda = \kappa \mu}{lambda=kappa* mu}.

  In this algorithm, the Method of Minimum Contrast is first used to find
  optimal values of the parameters \eqn{\kappa}{kappa}
  and \eqn{\sigma^2}{sigma^2}. Then the remaining parameter
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

  The Thomas process can be simulated, using \code{\link{rThomas}}.

  Homogeneous or inhomogeneous Thomas process models can also
  be fitted using the function \code{\link{kppm}}.

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
  Diggle, P. J., Besag, J. and Gleaves, J. T. (1976)
  Statistical analysis of spatial point patterns by
  means of distance methods. \emph{Biometrics} \bold{32} 659--667.

  Moller, J. and Waagepetersen, R. (2003).
  Statistical Inference and Simulation for Spatial Point Processes.
  Chapman and Hall/CRC, Boca Raton.

  Thomas, M. (1949) A generalisation of Poisson's binomial limit for use
  in ecology. \emph{Biometrika} \bold{36}, 18--25.

  Waagepetersen, R. (2007)
  An estimating function approach to inference for
  inhomogeneous Neyman-Scott processes.
  \emph{Biometrics} \bold{63}, 252--258.
}
\author{
  Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
}
\seealso{
  \code{\link{thomas.estK}}
  \code{\link{mincontrast}},
  \code{\link{pcf}},
  \code{\link{rThomas}} to simulate the fitted model.
}
\examples{
    data(redwood)
    u <- thomas.estpcf(redwood, c(kappa=10, sigma2=0.1))
    u
    plot(u, legendpos="topright")
    u2 <- thomas.estpcf(redwood, c(kappa=10, sigma2=0.1),
          pcfargs=list(stoyan=0.12))
}
\keyword{spatial}
\keyword{models}
