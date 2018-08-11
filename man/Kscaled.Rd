\name{Kscaled}
\alias{Kscaled}
\alias{Lscaled}
\title{Locally Scaled K-function}
\description{
  Estimates the locally-rescaled \eqn{K}-function of
  a point process.
}
\usage{
  Kscaled(X, lambda=NULL, \dots, r = NULL, breaks = NULL,
    rmax = 2.5, 
    correction=c("border", "isotropic", "translate"),
    renormalise=FALSE, normpower=1,
    sigma=NULL, varcov=NULL)

  Lscaled(\dots)
}
\arguments{
  \item{X}{
    The observed data point pattern,
    from which an estimate of the locally scaled \eqn{K} function
    will be computed.
    An object of class \code{"ppp"}
    or in a format recognised by \code{\link{as.ppp}()}.
  }
  \item{lambda}{
    Optional.
    Values of the estimated intensity function.
    Either a vector giving the intensity values
    at the points of the pattern \code{X},
    a pixel image (object of class \code{"im"}) giving the
    intensity values at all locations, a \code{function(x,y)} which
    can be evaluated to give the intensity value at any location,
    or a fitted point process model (object of class \code{"ppm"}).
  }
  \item{\dots}{
    Arguments passed from \code{Lscaled} to \code{Kscaled}
    and from \code{Kscaled} to \code{\link{density.ppp}}
    if \code{lambda} is omitted.
  }
  \item{r}{
    vector of values for the argument \eqn{r} at which
    the locally scaled \eqn{K} function
    should be evaluated. (These are rescaled distances.)
    Not normally given by the user; there is a sensible default.
  }
  \item{breaks}{
    This argument is for internal use only.
  }
  \item{rmax}{
    maximum value of the argument \eqn{r} that should be used.
    (This is the rescaled distance).
  }
  \item{correction}{
    A character vector containing any selection of the
    options \code{"border"}, 
    \code{"isotropic"}, \code{"Ripley"},
    \code{"translate"}, \code{"translation"},
    \code{"none"} or \code{"best"}.
    It specifies the edge correction(s) to be applied.
    Alternatively \code{correction="all"} selects all options.
  }
  \item{renormalise}{
    Logical. Whether to renormalise the estimate. See Details.
  }
  \item{normpower}{
    Integer (usually either 1 or 2).
    Normalisation power. See Details.
  }
  \item{sigma,varcov}{
    Optional arguments passed to  \code{\link{density.ppp}}
    to control the smoothing bandwidth, when \code{lambda} is
    estimated by kernel smoothing.
  }
}
\value{
  An object of class \code{"fv"} (see \code{\link{fv.object}}).
  
  Essentially a data frame containing at least the following columns,
  \item{r}{the vector of values of the argument \eqn{r} 
    at which the pair correlation function \eqn{g(r)} has been  estimated
  }
  \item{theo}{vector of values of \eqn{\pi r^2}{pi * r^2},
    the theoretical value of \eqn{K_{\rm scaled}(r)}{Kscaled(r)}
    for an inhomogeneous Poisson process
  }
  and containing additional columns
  according to the choice specified in the \code{correction}
  argument. The additional columns are named
  \code{border}, \code{trans} and \code{iso}
  and give the estimated values of 
  \eqn{K_{\rm scaled}(r)}{Kscaled(r)}
  using the border correction, translation correction,
  and Ripley isotropic correction, respectively.
}
\details{
  \code{Kscaled} computes an estimate of the \eqn{K} function
  for a locally scaled point process.
  \code{Lscaled} computes the corresponding \eqn{L} function
  \eqn{L(r) = \sqrt{K(r)/\pi}}{L(r) = sqrt(K(r)/pi)}.

  Locally scaled point processes are a class of models
  for inhomogeneous point patterns, introduced by Hahn et al (2003).
  They include inhomogeneous Poisson processes, and many other models.

  The template \eqn{K} function of a locally-scaled process is a counterpart
  of the ``ordinary'' Ripley \eqn{K} function, in which
  the distances between points of the process are measured
  on a spatially-varying scale (such that the locally rescaled
  process has unit intensity).

  The template \eqn{K} function is an indicator of interaction
  between the points. For an inhomogeneous Poisson process, the
  theoretical template \eqn{K} function is approximately equal
  to \eqn{K(r) = \pi r^2}{K(r) = pi * r^2}.
  Values \eqn{K_{\rm scaled}(r) > \pi r^2}{Kscaled(r) > pi * r^2}
  are suggestive of clustering.

  \code{Kscaled} computes an estimate of the template \eqn{K} function
  and \code{Lscaled} computes the corresponding \eqn{L} function
  \eqn{L(r) = \sqrt{K(r)/\pi}}{L(r) = sqrt(K(r)/pi)}.
  
  The locally scaled interpoint distances are computed
  using an approximation proposed by Hahn (2007). The Euclidean
  distance between two points is multiplied by the average of the
  square roots of the intensity values at the two points.

  The argument \code{lambda} should supply the
  (estimated) values of the intensity function \eqn{\lambda}{lambda}.
  It may be either
  \describe{
    \item{a numeric vector}{
      containing the values
      of the intensity function at the points of the pattern \code{X}.
    }
    \item{a pixel image}{
      (object of class \code{"im"})
      assumed to contain the values of the intensity function
      at all locations in the window. 
    }
    \item{a function}{
      which can be evaluated to give values of the intensity at
      any locations.
    }
    \item{omitted:}{
      if \code{lambda} is omitted, then it will be estimated using
      a `leave-one-out' kernel smoother.
    }
  }
  If \code{lambda} is a numeric vector, then its length should
  be equal to the number of points in the pattern \code{X}.
  The value \code{lambda[i]} is assumed to be the 
  the (estimated) value of the intensity
  \eqn{\lambda(x_i)}{lambda(x[i])} for
  the point \eqn{x_i}{x[i]} of the pattern \eqn{X}.
  Each value must be a positive number; \code{NA}'s are not allowed.

  If \code{lambda} is a pixel image, the domain of the image should
  cover the entire window of the point pattern. If it does not (which
  may occur near the boundary because of discretisation error),
  then the missing pixel values 
  will be obtained by applying a Gaussian blur to \code{lambda} using
  \code{\link{blur}}, then looking up the values of this blurred image
  for the missing locations. 
  (A warning will be issued in this case.)

  If \code{lambda} is a function, then it will be evaluated in the
  form \code{lambda(x,y)} where \code{x} and \code{y} are vectors
  of coordinates of the points of \code{X}. It should return a numeric
  vector with length equal to the number of points in \code{X}.

  If \code{lambda} is omitted, then it will be estimated using
  a `leave-one-out' kernel smoother,
  as described in Baddeley, \Moller
  and Waagepetersen (2000).  The estimate \code{lambda[i]} for the
  point \code{X[i]} is computed by removing \code{X[i]} from the
  point pattern, applying kernel smoothing to the remaining points using
  \code{\link{density.ppp}}, and evaluating the smoothed intensity
  at the point \code{X[i]}. The smoothing kernel bandwidth is controlled
  by the arguments \code{sigma} and \code{varcov}, which are passed to
  \code{\link{density.ppp}} along with any extra arguments.
  
  If \code{renormalise=TRUE}, the estimated intensity \code{lambda}
  is multiplied by \eqn{c^(normpower/2)} before performing other calculations,
  where \eqn{c = area(W)/sum[i] (1/lambda(x[i]))}. This 
  renormalisation has about the same effect as in \code{\link{Kinhom}},
  reducing the variability and bias of the estimate
  in small samples and in cases of very strong inhomogeneity.

  Edge corrections are used to correct bias in the estimation
  of \eqn{K_{\rm scaled}}{Kscaled}. First the interpoint distances are 
  rescaled, and then edge corrections are applied as in \code{\link{Kest}}.
  See \code{\link{Kest}} for details of the edge corrections
  and the options for the argument \code{correction}.
  
  The pair correlation function can also be applied to the
  result of \code{Kscaled}; see \code{\link{pcf}} and \code{\link{pcf.fv}}.
}
\references{
  Baddeley, A.,
  \Moller, J. and Waagepetersen, R. (2000)
  Non- and semiparametric estimation of interaction in
  inhomogeneous point patterns.
  \emph{Statistica Neerlandica} \bold{54}, 329--350.

  Hahn, U. (2007)
  \emph{Global and Local Scaling in the
    Statistics of Spatial Point Processes}. Habilitationsschrift,
  Universitaet Augsburg.
  
  Hahn, U., Jensen, E.B.V., van Lieshout, M.N.M. and Nielsen, L.S. (2003)
  Inhomogeneous spatial point processes by location-dependent scaling.
  \emph{Advances in Applied Probability} \bold{35}, 319--336.

  \Prokesova, M.,
  Hahn, U. and Vedel Jensen, E.B. (2006)
  Statistics for locally scaled point patterns.
  In A. Baddeley, P. Gregori, J. Mateu, R. Stoica and D. Stoyan (eds.)
  \emph{Case Studies in Spatial Point Pattern Modelling}.
  Lecture Notes in Statistics 185. New York: Springer Verlag.
  Pages 99--123.
}
\seealso{
  \code{\link{Kest}},
  \code{\link{pcf}}
}
\examples{
  data(bronzefilter)
  X <- unmark(bronzefilter)
  K <- Kscaled(X)
  fit <- ppm(X, ~x)
  lam <- predict(fit)
  K <- Kscaled(X, lam)
}
\author{Ute Hahn,
  \adrian
  
  
  and \rolf
  
}
\keyword{spatial}
\keyword{nonparametric}
