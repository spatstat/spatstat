\name{clusterfit}
\alias{clusterfit}
\title{Fit Cluster or Cox Point Process Model via Minimum Contrast}
\description{
  Fit a homogeneous or inhomogeneous cluster process or
  Cox point process model to a point pattern by the Method of Minimum Contrast.
}
\usage{
clusterfit(X, clusters, lambda = NULL, startpar = NULL, \dots,
           q = 1/4, p = 2, rmin = NULL, rmax = NULL,
           ctrl=list(q=q, p=p, rmin=rmin, rmax=rmax),
           statistic = NULL, statargs = NULL, algorithm="Nelder-Mead")
}
\arguments{
  \item{X}{
    Data to which the cluster or Cox model will be fitted.
    Either a point pattern or a summary statistic.
    See Details.
  }
  \item{clusters}{
    Character string determining the cluster or Cox model.
    Partially matched.
    Options are \code{"Thomas"}, \code{"MatClust"},
    \code{"Cauchy"}, \code{"VarGamma"} and \code{"LGCP"}.
  }
  \item{lambda}{
    Optional. An estimate of the intensity of the point process.
    Either a single numeric specifying a constant intensity,
    a pixel image (object of class \code{"im"}) giving the
    intensity values at all locations, a fitted point process model
    (object of class \code{"ppm"} or \code{"kppm"})
    or a \code{function(x,y)} which
    can be evaluated to give the intensity value at any location.
  }
  \item{startpar}{
    Vector of initial values of the parameters of the
    point process mode. If \code{X} is a point pattern sensible defaults
    are used. Otherwise rather arbitrary values are used.
  }
  \item{q,p}{
    Optional. Exponents for the contrast criterion.
    See \code{\link{mincontrast}}.
  }
  \item{rmin, rmax}{
    Optional. The interval of \eqn{r} values for the contrast criterion.
    See \code{\link{mincontrast}}.
  }
  \item{ctrl}{
    Optional. Named list containing values of the parameters
    \code{q,p,rmin,rmax}.
  }
  \item{\dots}{
    Additional arguments passed to \code{\link{mincontrast}.}
  }
  \item{statistic}{
    Optional. Name of the summary statistic to be used
    for minimum contrast estimation: either \code{"K"} or \code{"pcf"}.
  }
  \item{statargs}{
    Optional list of arguments to be used when calculating
    the \code{statistic}. See Details.
  }
  \item{algorithm}{
    Character string determining the mathematical optimisation algorithm
    to be used by \code{\link[stats]{optim}}. See
    the argument \code{method} of \code{\link[stats]{optim}}.
  }
}
\details{
  This function fits the clustering parameters of a cluster or Cox point
  process model by the Method of Minimum Contrast, that is, by
  matching the theoretical \eqn{K}-function of the model to the
  empirical \eqn{K}-function of the data, as explained in
  \code{\link{mincontrast}}.
  
  If \code{statistic="pcf"} (or \code{X} appears to be an
  estimated pair correlation function) then instead of using the
  \eqn{K}-function, the algorithm will use the pair correlation
  function.

  If \code{X} is a point pattern of class \code{"ppp"} an estimate of
  the summary statistic specfied by \code{statistic} (defaults to
  \code{"K"}) is first computed before minimum contrast estimation is
  carried out as described above. In this case the argument
  \code{statargs} can be used for controlling the summary statistic
  estimation. The precise algorithm for computing the summary statistic
  depends on whether the intensity specification (\code{lambda}) is:

  \describe{
    \item{homogeneous:}{
      If \code{lambda} is \code{NUll} or a single numeric the pattern is
      considered homogeneous and either \code{\link{Kest}} or
      \code{\link{pcf}} is invoked. In this case \code{lambda} is
      \bold{not} used for anything when estimating the summary statistic.
    }
    \item{inhomogeneous:}{

      If \code{lambda} is a pixel image (object of class \code{"im"}),
      a fitted point process model (object of class \code{"ppm"} or
      \code{"kppm"}) or a \code{function(x,y)} the pattern is considered
      inhomogeneous. In this case either \code{\link{Kinhom}} or
      \code{\link{pcfinhom}} is invoked with \code{lambda} as an
      argument.
      
    }
  }

  After the clustering parameters of the model have been estimated by
  minimum contrast \code{lambda} (if non-null) is used to compute the
  additional model parameter \eqn{\mu}{\mu}.

  The algorithm parameters \code{q,p,rmax,rmin} are described in the
  help for \code{\link{mincontrast}}. They may be provided either
  as individually-named arguments, or as entries in the list
  \code{ctrl}. The individually-named arguments \code{q,p,rmax,rmin}
  override the entries in the list \code{ctrl}.
}
\value{
  An object of class \code{"minconfit"}. There are methods for printing
  and plotting this object. See \code{\link{mincontrast}}.
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
\author{
  \spatstatAuthors.
}
\seealso{
  \code{\link{kppm}}
}
\examples{
  fit <- clusterfit(redwood, "Thomas")
  fit
  if(interactive()){
    plot(fit)
  }
  K <- Kest(redwood)
  fit2 <- clusterfit(K, "MatClust")
}
\keyword{spatial}
\keyword{models}
