\name{rhohat}  
\alias{rhohat}
\alias{rhohat.ppp}
\alias{rhohat.quad}
\alias{rhohat.ppm}
\alias{rhohat.lpp}
\alias{rhohat.lppm}
\concept{Resource Selection Function}
\concept{Prospectivity}
\title{
  Smoothing Estimate of Intensity as Function of a Covariate
}
\description{
  Computes a smoothing estimate of the intensity of a point process,
  as a function of a (continuous) spatial covariate.
}
\usage{
rhohat(object, covariate, ...)

\method{rhohat}{ppp}(object, covariate, ...,
       baseline=NULL, weights=NULL,
       method=c("ratio", "reweight", "transform"),
       horvitz=FALSE,
       smoother=c("kernel", "local"),
       subset=NULL, 
       dimyx=NULL, eps=NULL,
       n = 512, bw = "nrd0", adjust=1, from = NULL, to = NULL,
       bwref=bw,
       covname, confidence=0.95, positiveCI)

\method{rhohat}{quad}(object, covariate, ...,
       baseline=NULL, weights=NULL,
       method=c("ratio", "reweight", "transform"),
       horvitz=FALSE,
       smoother=c("kernel", "local"),
       subset=NULL, 
       dimyx=NULL, eps=NULL,
       n = 512, bw = "nrd0", adjust=1, from = NULL, to = NULL,
       bwref=bw,
       covname, confidence=0.95, positiveCI)

\method{rhohat}{ppm}(object, covariate, ...,
       weights=NULL,
       method=c("ratio", "reweight", "transform"),
       horvitz=FALSE,
       smoother=c("kernel", "local"),
       subset=NULL, 
       dimyx=NULL, eps=NULL,
       n = 512, bw = "nrd0", adjust=1, from = NULL, to = NULL,
       bwref=bw,
       covname, confidence=0.95, positiveCI)

\method{rhohat}{lpp}(object, covariate, ...,
       weights=NULL,
       method=c("ratio", "reweight", "transform"),
       horvitz=FALSE,
       smoother=c("kernel", "local"),
       subset=NULL, 
       nd=1000, eps=NULL, random=TRUE, 
       n = 512, bw = "nrd0", adjust=1, from = NULL, to = NULL,
       bwref=bw,
       covname, confidence=0.95, positiveCI)

\method{rhohat}{lppm}(object, covariate, ...,
       weights=NULL,
       method=c("ratio", "reweight", "transform"),
       horvitz=FALSE,
       smoother=c("kernel", "local"),
       subset=NULL, 
       nd=1000, eps=NULL, random=TRUE, 
       n = 512, bw = "nrd0", adjust=1, from = NULL, to = NULL,
       bwref=bw,
       covname, confidence=0.95, positiveCI)
}
\arguments{
  \item{object}{
    A point pattern (object of class \code{"ppp"} or \code{"lpp"}),
    a quadrature scheme (object of class \code{"quad"})
    or a fitted point process model (object of class \code{"ppm"}
    or \code{"lppm"}).
  }
  \item{covariate}{
    Either a \code{function(x,y)} or a pixel image (object of
    class \code{"im"}) providing the values of the covariate at any
    location.
    Alternatively one of the strings \code{"x"} or \code{"y"}
    signifying the Cartesian coordinates.
  }
  \item{weights}{
    Optional weights attached to the data points.
    Either a numeric vector of weights for each data point,
    or a pixel image (object of class \code{"im"}) or
    a \code{function(x,y)} providing the weights.
  }
  \item{baseline}{
    Optional baseline for intensity function.
    A \code{function(x,y)} or a pixel image (object of
    class \code{"im"}) providing the values of the baseline at any
    location. 
  }
  \item{method}{
    Character string determining the smoothing method. See Details.
  }
  \item{horvitz}{
    Logical value indicating whether to use Horvitz-Thompson weights.
    See Details.
  }
  \item{smoother}{
    Character string determining the smoothing algorithm. See Details.
  }
  \item{subset}{
    Optional. A spatial window (object of class \code{"owin"})
    specifying a subset of the data, from which the estimate should
    be calculated.
  }
  \item{dimyx,eps,nd,random}{
    Arguments controlling the pixel
    resolution at which the covariate will be evaluated.
    See Details.
  }
  \item{bw}{
    Smoothing bandwidth or bandwidth rule
    (passed to \code{\link{density.default}}).
  }
  \item{adjust}{
    Smoothing bandwidth adjustment factor
    (passed to \code{\link{density.default}}).
  }
  \item{n, from, to}{
    Arguments passed to \code{\link{density.default}} to
    control the number and range of values at which the function
    will be estimated.
  }
  \item{bwref}{
    Optional. An alternative value of \code{bw} to use when smoothing
    the reference density (the density of the covariate values
    observed at all locations in the window).
  }
  \item{\dots}{
    Additional arguments passed to \code{\link{density.default}}
    or \code{\link[locfit]{locfit}}.
  }
  \item{covname}{
    Optional. Character string to use as the name of the covariate.
  }
  \item{confidence}{
    Confidence level for confidence intervals.
    A number between 0 and 1.
  }
  \item{positiveCI}{
    Logical value.
    If \code{TRUE}, confidence limits are always positive numbers;
    if \code{FALSE}, the lower limit of the
    confidence interval may sometimes be negative.
    Default is \code{FALSE} if \code{smoother="kernel"}
    and \code{TRUE} if \code{smoother="local"}.
    See Details.
  }
}
\details{
  This command estimates the relationship between
  point process intensity and a given spatial covariate.
  Such a relationship is sometimes called a
  \emph{resource selection function} (if the points are organisms
  and the covariate is a descriptor of habitat) or
  a \emph{prospectivity index} (if the points are mineral deposits
  and the covariate is a geological variable). 
  This command uses a nonparametric smoothing method which does not assume a
  particular form for the relationship.  
  
  If \code{object} is a point pattern, and \code{baseline} is missing or
  null, this command assumes that \code{object} is a realisation of a
  Poisson point process with intensity function
  \eqn{\lambda(u)}{lambda(u)} of the form
  \deqn{\lambda(u) = \rho(Z(u))}{lambda(u) = rho(Z(u))}
  where \eqn{Z} is the spatial
  covariate function given by \code{covariate}, and
  \eqn{\rho(z)}{rho(z)} is a function to be estimated.  This command
  computes estimators of \eqn{\rho(z)}{rho(z)} proposed by Baddeley and
  Turner (2005) and Baddeley et al (2012).

  The covariate \eqn{Z} must have continuous values.

  If \code{object} is a point pattern, and \code{baseline} is given,
  then the intensity function is assumed to be
  \deqn{\lambda(u) = \rho(Z(u)) B(u)}{lambda(u) = rho(Z(u)) * B(u)}
  where \eqn{B(u)} is the baseline intensity at location \eqn{u}.
  A smoothing estimator of the relative intensity  \eqn{\rho(z)}{rho(z)}
  is computed.

  If \code{object} is a fitted point process model, suppose \code{X} is
  the original data point pattern to which the model was fitted. Then
  this command assumes \code{X} is a realisation of a Poisson point
  process with intensity function of the form
  \deqn{
    \lambda(u) = \rho(Z(u)) \kappa(u)
  }{
    lambda(u) = rho(Z(u)) * kappa(u)
  }
  where \eqn{\kappa(u)}{kappa(u)} is the intensity of the fitted model
  \code{object}. A smoothing estimator of \eqn{\rho(z)}{rho(z)} is computed.

  The estimation procedure is determined by the character strings
  \code{method} and \code{smoother} and the argument \code{horvitz}.
  The estimation procedure involves computing several density estimates
  and combining them. 
  The algorithm used to compute density estimates is 
  determined by \code{smoother}:
  \itemize{
    \item If \code{smoother="kernel"},
    each the smoothing procedure is based on
    fixed-bandwidth kernel density estimation,
    performed by \code{\link{density.default}}.
    \item If \code{smoother="local"}, the smoothing procedure
    is based on local likelihood density estimation, performed by
    \code{\link[locfit]{locfit}}.
  }
  The \code{method} determines how the density estimates will be
  combined to obtain an estimate of \eqn{\rho(z)}{rho(z)}:
  \itemize{
    \item
    If \code{method="ratio"}, then \eqn{\rho(z)}{rho(z)} is
    estimated by the ratio of two density estimates.
    The numerator is a (rescaled) density estimate obtained by
    smoothing the values \eqn{Z(y_i)}{Z(y[i])} of the covariate
    \eqn{Z} observed at the data points \eqn{y_i}{y[i]}. The denominator
    is a density estimate of the reference distribution of \eqn{Z}.
    \item
    If \code{method="reweight"}, then \eqn{\rho(z)}{rho(z)} is
    estimated by applying density estimation to the 
    values \eqn{Z(y_i)}{Z(y[i])} of the covariate
    \eqn{Z} observed at the data points \eqn{y_i}{y[i]},
    with weights inversely proportional to the reference density of
    \eqn{Z}.
    \item 
    If \code{method="transform"},
    the smoothing method is variable-bandwidth kernel
    smoothing, implemented by applying the Probability Integral Transform
    to the covariate values, yielding values in the range 0 to 1,
    then applying edge-corrected density estimation on the interval
    \eqn{[0,1]}, and back-transforming.
  }
  If \code{horvitz=TRUE}, then the calculations described above
  are modified by using Horvitz-Thompson weighting.
  The contribution to the numerator from 
  each data point is weighted by the reciprocal of the
  baseline value or fitted intensity value at that data point;
  and a corresponding adjustment is made to the denominator.
  
  The covariate will be evaluated on a fine grid of locations,
  with spatial resolution controlled by the arguments
  \code{dimyx,eps,nd,random}. 
  In two dimensions (i.e.
  if \code{object} is of class \code{"ppp"}, \code{"ppm"} or
  \code{"quad"}) the arguments \code{dimyx, eps} are
  passed to \code{\link{as.mask}} to control the pixel
  resolution. On a linear network (i.e. if \code{object} is of class
  \code{"lpp"}) the argument \code{nd} specifies the
  total number of test locations on the linear
  network, \code{eps} specifies the linear separation between test
  locations, and \code{random} specifies whether the test locations
  have a randomised starting position.

  If the argument \code{weights} is present, then the contribution
  from each data point \code{X[i]} to the estimate of \eqn{\rho}{rho} is
  multiplied by \code{weights[i]}.

  If the argument \code{subset} is present, then the calculations are
  performed using only the data inside this spatial region.
    
  Pointwise confidence intervals for the true value of \eqn{\rho(z)}
  are also calculated for each \eqn{z} using the central limit theorem,
  and will be plotted as grey shading.
  If \code{positiveCI=FALSE}, the lower limit of the confidence
  interval may sometimes be negative, because the confidence intervals
  are based on a normal approximation to the estimate of \eqn{\rho(z)}.
  If \code{positiveCI=TRUE}, the confidence limits are always
  positive, because the confidence interval is based on a normal
  approximation to the estimate of \eqn{\log(\rho(z))}{log(\rho(z))}.
  For consistency with earlier versions, the default is
  \code{positiveCI=FALSE} for \code{smoother="kernel"}
  and \code{positiveCI=TRUE} for \code{smoother="local"}.
}
\value{
  A function value table (object of class \code{"fv"})
  containing the estimated values of \eqn{\rho}{rho} for a sequence
  of values of \eqn{Z}.
  Also belongs to the class \code{"rhohat"}
  which has special methods for \code{print}, \code{plot}
  and \code{predict}.
}
\section{Categorical and discrete covariates}{
  This technique assumes that the covariate has continuous values.
  It is not applicable to covariates with categorical (factor) values
  or discrete values such as small integers.
  For a categorical covariate, use
  \code{\link{intensity.quadratcount}} applied to the result of
  \code{\link{quadratcount}(X, tess=covariate)}.
}
\references{
  Baddeley, A., Chang, Y.-M., Song, Y. and Turner, R. (2012)
  Nonparametric estimation of the dependence of a point
  process on spatial covariates.
  \emph{Statistics and Its Interface} \bold{5} (2), 221--236.
  
  Baddeley, A. and Turner, R. (2005)
  Modelling spatial point patterns in R.
  In: A. Baddeley, P. Gregori, J. Mateu, R. Stoica, and D. Stoyan,
  editors, \emph{Case Studies in Spatial Point Pattern Modelling},
  Lecture Notes in Statistics number 185. Pages 23--74.
  Springer-Verlag, New York, 2006. 
  ISBN: 0-387-28311-0.  
}
\author{
  \adrian, 
  Ya-Mei Chang, Yong Song, 
  and \rolf.
}
\seealso{
  \code{\link{rho2hat}},
  \code{\link{methods.rhohat}},
  \code{\link{parres}}.

  See \code{\link{ppm}} for a parametric method for the same problem.
}
\examples{
  X <-  rpoispp(function(x,y){exp(3+3*x)})
  rho <- rhohat(X, "x")
  rho <- rhohat(X, function(x,y){x})
  plot(rho)
  curve(exp(3+3*x), lty=3, col=2, add=TRUE)

  rhoB <- rhohat(X, "x", method="reweight")
  rhoC <- rhohat(X, "x", method="transform")

  \testonly{rh <- rhohat(X, "x", dimyx=32)}

  fit <- ppm(X, ~x)
  rr <- rhohat(fit, "y")

# linear network
  Y <- runiflpp(30, simplenet)
  rhoY <- rhohat(Y, "y")
}
\keyword{spatial}
\keyword{models}

