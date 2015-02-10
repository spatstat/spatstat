\name{clusterfit}
\alias{clusterfit}
\title{Fit Cluster or Cox Point Process Model via Minimum Contrast}
\description{
  Fit a homogeneous or inhomogeneous cluster process or
  Cox point process model to a point pattern by the Method of Minimum Contrast.
}
\usage{
clusterfit(X, clusters, lambda = NULL, startpar = NULL,
           q = 1/4, p = 2, rmin = NULL, rmax = NULL, \dots,
           statistic = NULL, statargs = NULL)
}
\arguments{
  \item{X}{
    Data to which the cluster or Cox model will be fitted.
    Either a point pattern or a summary statistic.
    See Details.
  }
  \item{clusters}{
    Character string determining the cluster model.
    Partially matched.
    Options are \code{"Thomas"}, \code{"MatClust"},
    \code{"Cauchy"}, \code{"VarGamma"} and \code{"LGCP"}.
  }
  \item{lambda}{
    Optional. An estimate of the intensity of the point process.
  }
  \item{startpar}{
    Vector of initial values of the parameters of the
    point process mode.
  }
  \item{q,p}{
    Optional. Exponents for the contrast criterion.
  }
  \item{rmin, rmax}{
    Optional. The interval of \eqn{r} values for the contrast criterion.
  }
  \item{\dots}{
    Additional arguments passed to \code{\link{mincontrast}.}
  }
  \item{statistic}{
    Name of the summary statistic to be used
    for minimum contrast estimation: either \code{"K"} or \code{"pcf"}.
  }
  \item{statargs}{
    Optional list of arguments to be used when calculating
    the \code{statistic}. See Details.
  }
}
\details{
  This function fits the clustering parameters of a cluster or Cox point
  process model by the Method of Minimum Contrast, that is, by
  matching the theoretical \eqn{K}-function of the model to the
  empirical \eqn{K}-function of the data, as explained in
  \code{\link{mincontrast}}....
}
\value{
  An object of class \code{"minconfit"}. There are methods for printing
  and plotting this object. It contains the following components:
  \item{par }{Vector of fitted parameter values.}
  \item{fit }{Function value table (object of class \code{"fv"})
    containing the observed values of the summary statistic
    (\code{observed}) and the theoretical values of the summary
    statistic computed from the fitted model parameters.
  }
  \item{opt }{The return value from the optimizer \code{\link{optim}}.}
  \item{crtl }{The control parameters of the algorithm.}
  \item{info }{List of explanatory strings.}
}
\references{
  Diggle, P.J. and Gratton, R.J. (1984)
  Monte Carlo methods of inference for implicit statistical models.
  \emph{Journal of the Royal Statistical Society, series B}
  \bold{46}, 193 -- 212.

  \ifelse{latex}{\out{M\o ller}}{Moller}, J. and Waagepetersen, R. (2003).
  Statistical Inference and Simulation for Spatial Point Processes.
  Chapman and Hall/CRC, Boca Raton.

  Waagepetersen, R. (2007).
  An estimating function approach to inference for
  inhomogeneous Neyman-Scott processes.
  \emph{Biometrics} \bold{63} (2007) 252--258.
}  
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/},
  Rolf Turner
  \email{r.turner@auckland.ac.nz}
  and Ege Rubak
  \email{rubak@math.aau.dk}
}
\seealso{
  \code{\link{kppm}}
}
\keyword{spatial}
\keyword{models}
