\name{Kinhom}
\alias{Kinhom}
\title{Inhomogeneous K-function}
\description{
  Estimates the inhomogeneous \eqn{K} function of
  a non-stationary point pattern.
}
\usage{
  Kinhom(X, lambda=NULL, \dots, r = NULL, breaks = NULL,
    correction=c("border", "bord.modif", "isotropic", "translate"),
    renormalise=TRUE,
    normpower=1,
    update=TRUE,
    leaveoneout=TRUE,
    nlarge = 1000,
    lambda2=NULL, reciplambda=NULL, reciplambda2=NULL,
    sigma=NULL, varcov=NULL)
}
\arguments{
  \item{X}{
    The observed data point pattern,
    from which an estimate of the inhomogeneous \eqn{K} function
    will be computed.
    An object of class \code{"ppp"}
    or in a format recognised by \code{\link{as.ppp}()}
  }
  \item{lambda}{
    Optional.
    Values of the estimated intensity function.
    Either a vector giving the intensity values
    at the points of the pattern \code{X},
    a pixel image (object of class \code{"im"}) giving the
    intensity values at all locations, a fitted point process model
    (object of class \code{"ppm"} or \code{"kppm"})
    or a \code{function(x,y)} which
    can be evaluated to give the intensity value at any location.
  }
  \item{\dots}{
    Extra arguments. Ignored if \code{lambda} is present.
    Passed to \code{\link{density.ppp}} if \code{lambda} is omitted.
  }
  \item{r}{
    vector of values for the argument \eqn{r} at which
    the inhomogeneous \eqn{K} function
    should be evaluated. Not normally given by the user;
    there is a sensible default.
  }
  \item{breaks}{
    This argument is for internal use only.
  }
  \item{correction}{
    A character vector containing any selection of the
    options \code{"border"}, \code{"bord.modif"},
    \code{"isotropic"}, \code{"Ripley"}, \code{"translate"},
    \code{"translation"},
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
  \item{update}{
    Logical. If \code{lambda} is a fitted model
    (class \code{"ppm"}, \code{"kppm"} or \code{"dppm"})
    and \code{update=TRUE} (the default),
    the model will first be refitted to the data \code{X}
    (using \code{\link{update.ppm}} or \code{\link{update.kppm}})
    before the fitted intensity is computed.
    If \code{update=FALSE}, the fitted intensity of the
    model will be computed without re-fitting it to \code{X}.
  }
  \item{leaveoneout}{
    Logical value (passed to \code{\link{density.ppp}} or
    \code{\link{fitted.ppm}}) specifying whether to use a
    leave-one-out rule when calculating the intensity.
  }
  \item{nlarge}{
    Optional. Efficiency threshold.
    If the number of points exceeds \code{nlarge}, then only the
    border correction will be computed, using a fast algorithm.
  }
  \item{lambda2}{
    Advanced use only.
    Matrix containing estimates of the products
    \eqn{\lambda(x_i)\lambda(x_j)}{lambda(x[i]) * lambda(x[j])}
    of the intensities at each pair of data points 
    \eqn{x_i}{x[i]} and \eqn{x_j}{x[j]}. 
  }
  \item{reciplambda}{
    Alternative to \code{lambda}.
    Values of the estimated \emph{reciprocal} \eqn{1/\lambda}{1/lambda}
    of the intensity function.
    Either a vector giving the reciprocal intensity values
    at the points of the pattern \code{X},
    a pixel image (object of class \code{"im"}) giving the
    reciprocal intensity values at all locations,
    or a \code{function(x,y)} which can be evaluated to give the
    reciprocal intensity value at any location.
  }
  \item{reciplambda2}{
    Advanced use only. Alternative to \code{lambda2}.
    A matrix giving values of the estimated \emph{reciprocal products}
    \eqn{1/\lambda(x_i)\lambda(x_j)}{1/(lambda(x[i]) * lambda(x[j]))}
    of the intensities at each pair of data points 
    \eqn{x_i}{x[i]} and \eqn{x_j}{x[j]}. 
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
    at which \eqn{K_{\rm inhom}(r)}{Kinhom(r)} has been estimated
  }
  \item{theo}{vector of values of \eqn{\pi r^2}{pi * r^2},
    the theoretical value of \eqn{K_{\rm inhom}(r)}{Kinhom(r)}
    for an inhomogeneous Poisson process
  }
  and containing additional columns
  according to the choice specified in the \code{correction}
  argument. The additional columns are named
  \code{border}, \code{trans} and \code{iso}
  and give the estimated values of 
  \eqn{K_{\rm inhom}(r)}{Kinhom(r)}
  using the border correction, translation correction,
  and Ripley isotropic correction, respectively.
}
\details{
  This computes a generalisation of the \eqn{K} function
  for inhomogeneous point patterns, proposed by
  Baddeley, \ifelse{latex}{\out{M\o ller}}{Moller} and Waagepetersen (2000).
  
  The ``ordinary'' \eqn{K} function
  (variously known as the reduced second order moment function
  and Ripley's \eqn{K} function), is
  described under \code{\link{Kest}}. It is defined only
  for stationary point processes.
  
  The inhomogeneous \eqn{K} function
  \eqn{K_{\rm inhom}(r)}{Kinhom(r)}
  is a direct generalisation to nonstationary point processes.
  Suppose \eqn{x} is a point process with non-constant intensity
  \eqn{\lambda(u)}{lambda(u)} at each location \eqn{u}.
  Define \eqn{K_{\rm inhom}(r)}{Kinhom(r)} to be the expected
  value, given that \eqn{u} is a point of \eqn{x},
  of the sum of all terms
  \eqn{1/\lambda(x_j)}{1/lambda(x[j])}
  over all points \eqn{x_j}{x[j]}
  in the process separated from \eqn{u} by a distance less than \eqn{r}.
  This reduces to the ordinary \eqn{K} function if
  \eqn{\lambda()}{lambda()} is constant.
  If \eqn{x} is an inhomogeneous Poisson process with intensity
  function \eqn{\lambda(u)}{lambda(u)}, then
  \eqn{K_{\rm inhom}(r) = \pi r^2}{Kinhom(r) = pi * r^2}.
  
  Given a point pattern dataset, the 
  inhomogeneous \eqn{K} function can be estimated
  essentially by summing the values
  \eqn{1/(\lambda(x_i)\lambda(x_j))}{1/(lambda(x[i]) * lambda(x[j]))}
  for all pairs of points \eqn{x_i, x_j}{x[i], x[j]}
  separated by a distance less than \eqn{r}.

  This allows us to inspect a point pattern for evidence of 
  interpoint interactions after allowing for spatial inhomogeneity
  of the pattern. Values 
  \eqn{K_{\rm inhom}(r) > \pi r^2}{Kinhom(r) > pi * r^2}
  are suggestive of clustering.

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
    \item{a fitted point process model}{
      (object of class \code{"ppm"} or \code{"kppm"})
      whose fitted \emph{trend} can be used as the fitted intensity.
      (If \code{update=TRUE} the model will first be refitted to the
      data \code{X} before the trend is computed.)
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
  as described in Baddeley, \ifelse{latex}{\out{M\o ller}}{Moller}
  and Waagepetersen (2000).  The estimate \code{lambda[i]} for the
  point \code{X[i]} is computed by removing \code{X[i]} from the
  point pattern, applying kernel smoothing to the remaining points using
  \code{\link{density.ppp}}, and evaluating the smoothed intensity
  at the point \code{X[i]}. The smoothing kernel bandwidth is controlled
  by the arguments \code{sigma} and \code{varcov}, which are passed to
  \code{\link{density.ppp}} along with any extra arguments.
  
  Edge corrections are used to correct bias in the estimation
  of \eqn{K_{\rm inhom}}{Kinhom}.
  Each edge-corrected estimate of \eqn{K_{\rm inhom}(r)}{Kinhom(r)} is
  of the form
  \deqn{\widehat K_{\rm inhom}(r) = \sum_i \sum_j \frac{1\{d_{ij} \le
      r\} e(x_i,x_j,r)}{\lambda(x_i)\lambda(x_j)}
  }{K^inhom(r)= sum[i] sum[j] 1(d[i,j] <= r) * 
    e(x[i],x[j],r)/(lambda(x[i]) * lambda(x[j]))
  }
  where \eqn{d_{ij}}{d[i,j]} is the distance between points
  \eqn{x_i}{x[i]} and \eqn{x_j}{x[j]}, and
  \eqn{e(x_i,x_j,r)}{e(x[i],x[j],r)} is
  an edge correction factor. For the `border' correction,
  \deqn{
    e(x_i,x_j,r) =
    \frac{1(b_i > r)}{\sum_j 1(b_j > r)/\lambda(x_j)}
  }{
    1(b[i] > r)/(sum[j] 1(b[j] > r)/lambda(x[j]))
  }
  where \eqn{b_i}{b[i]} is the distance from \eqn{x_i}{x[i]}
  to the boundary of the window. For the `modified border'
  correction, 
  \deqn{
    e(x_i,x_j,r) =
    \frac{1(b_i > r)}{\mbox{area}(W \ominus r)}
  }{
    1(b[i] > r)/area(W [-] r)
  }
  where \eqn{W \ominus r}{W [-] r} is the eroded window obtained
  by trimming a margin of width \eqn{r} from the border of the original
  window.
  For the `translation' correction,
  \deqn{
    e(x_i,x_j,r) =
    \frac 1 {\mbox{area}(W \cap (W + (x_j - x_i)))}
  }{
    1/area(W intersect (W + x[j]-x[i]))
  }
  and for the `isotropic' correction,
  \deqn{
    e(x_i,x_j,r) =
    \frac 1 {\mbox{area}(W) g(x_i,x_j)}
  }{
    1/(area(W) g(x[i],x[j]))
  }
  where \eqn{g(x_i,x_j)}{g(x[i],x[j])} is the fraction of the
  circumference of the circle with centre \eqn{x_i}{x[i]} and radius
  \eqn{||x_i - x_j||}{||x[i]-x[j]||} which lies inside the window.
  
  If \code{renormalise=TRUE} (the default), then the estimates 
  are multiplied by \eqn{c^{\mbox{normpower}}}{c^normpower} where 
  \eqn{
    c = \mbox{area}(W)/\sum (1/\lambda(x_i)).
  }{
    c = area(W)/sum[i] (1/lambda(x[i])).
  }
  This rescaling reduces the variability and bias of the estimate
  in small samples and in cases of very strong inhomogeneity.
  The default value of \code{normpower} is 1 (for consistency with
  previous versions of \pkg{spatstat})
  but the most sensible value is 2, which would correspond to rescaling
  the \code{lambda} values so that
  \eqn{
    \sum (1/\lambda(x_i)) = \mbox{area}(W).
  }{
    sum[i] (1/lambda(x[i])) = area(W).
  }
  
  If the point pattern \code{X} contains more than about 1000 points,
  the isotropic and translation edge corrections can be computationally
  prohibitive. The computations for the border method are much faster,
  and are statistically efficient when there are large numbers of
  points. Accordingly, if the number of points in \code{X} exceeds
  the threshold \code{nlarge}, then only the border correction will be
  computed. Setting \code{nlarge=Inf} or \code{correction="best"}
  will prevent this from happening.
  Setting \code{nlarge=0} is equivalent to selecting only the border
  correction with \code{correction="border"}.

  The pair correlation function can also be applied to the
  result of \code{Kinhom}; see \code{\link{pcf}}.
}
\references{
  Baddeley, A.,
  \ifelse{latex}{\out{M\o ller}}{Moller}, J. and Waagepetersen, R. (2000)
  Non- and semiparametric estimation of interaction in
  inhomogeneous point patterns.
  \emph{Statistica Neerlandica} \bold{54}, 329--350.
}
\seealso{
  \code{\link{Kest}},
  \code{\link{pcf}}
}
\examples{
  # inhomogeneous pattern of maples
  X <- unmark(split(lansing)$maple)
  \testonly{
     sub <- sample(c(TRUE,FALSE), npoints(X), replace=TRUE, prob=c(0.1,0.9))
     X <- X[sub]
  }

  # (1) intensity function estimated by model-fitting
  # Fit spatial trend: polynomial in x and y coordinates
  fit <- ppm(X, ~ polynom(x,y,2), Poisson())
  # (a) predict intensity values at points themselves,
  #     obtaining a vector of lambda values
  lambda <- predict(fit, locations=X, type="trend")
  # inhomogeneous K function
  Ki <- Kinhom(X, lambda)
  plot(Ki)
  # (b) predict intensity at all locations,
  #     obtaining a pixel image
  lambda <- predict(fit, type="trend")
  Ki <- Kinhom(X, lambda)
  plot(Ki)

  # (2) intensity function estimated by heavy smoothing
  Ki <- Kinhom(X, sigma=0.1)
  plot(Ki)

  # (3) simulated data: known intensity function
  lamfun <- function(x,y) { 50 + 100 * x }
  # inhomogeneous Poisson process
  Y <- rpoispp(lamfun, 150, owin())
  # inhomogeneous K function
  Ki <- Kinhom(Y, lamfun)
  plot(Ki)

  # How to make simulation envelopes:
  #      Example shows method (2)
  \dontrun{
  smo <- density.ppp(X, sigma=0.1)
  Ken <- envelope(X, Kinhom, nsim=99,
                  simulate=expression(rpoispp(smo)),
                  sigma=0.1, correction="trans")
  plot(Ken)
  }
}
\author{Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  and Rolf Turner \email{r.turner@auckland.ac.nz}
  
}
\keyword{spatial}
\keyword{nonparametric}
