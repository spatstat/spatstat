\name{pcf.ppp}
\alias{pcf.ppp}
\title{Pair Correlation Function of Point Pattern}
\description{
  Estimates the pair correlation function of
  a point pattern using kernel methods.
}
\usage{
  \method{pcf}{ppp}(X, \dots, r = NULL, kernel="epanechnikov", bw=NULL,
                    stoyan=0.15,
                    correction=c("translate", "Ripley"),
                    divisor = c("r", "d"),
                    var.approx = FALSE,
                    domain=NULL,
                    ratio=FALSE)
}
\arguments{
  \item{X}{
    A point pattern (object of class \code{"ppp"}).
  }
  \item{r}{
    Vector of values for the argument \eqn{r} at which \eqn{g(r)} 
    should be evaluated. There is a sensible default.
  }
  \item{kernel}{
    Choice of smoothing kernel,
    passed to \code{density}.
  }
  \item{bw}{
    Bandwidth for smoothing kernel, passed to \code{density}.
  }
  \item{\dots}{
    Other arguments passed to the kernel density estimation 
    function \code{density}.
  }
  \item{stoyan}{
    Bandwidth coefficient; see Details.
  }
  \item{correction}{
    Choice of edge correction.
  }
  \item{divisor}{
    Choice of divisor in the estimation formula:
    either \code{"r"} (the default) or \code{"d"}. See Details.
  }
  \item{var.approx}{
    Logical value indicating whether to compute an analytic
    approximation to the variance of the estimated pair correlation.
  }
  \item{domain}{
    Optional. Calculations will be restricted to this subset
    of the window. See Details.
  }
  \item{ratio}{
    Logical. 
    If \code{TRUE}, the numerator and denominator of
    each edge-corrected estimate will also be saved,
    for use in analysing replicated point patterns.
  }
} 
\value{
  A function value table
  (object of class \code{"fv"}).
  Essentially a data frame containing the variables
  \item{r}{the vector of values of the argument \eqn{r} 
    at which the pair correlation function \eqn{g(r)} has been  estimated
  }
  \item{theo}{vector of values equal to 1,
    the theoretical value of \eqn{g(r)} for the Poisson process
  }
  \item{trans}{vector of values of \eqn{g(r)}
    estimated by translation correction
  }
  \item{iso}{vector of values of \eqn{g(r)}
    estimated by Ripley isotropic correction
  }
  \item{v}{vector of approximate values of the variance of
    the estimate of \eqn{g(r)}
  }
  as required.

  If \code{ratio=TRUE} then the return value also has two
  attributes called \code{"numerator"} and \code{"denominator"}
  which are \code{"fv"} objects
  containing the numerators and denominators of each
  estimate of \eqn{g(r)}.

  The return value also has an attribute \code{"bw"} giving the
  smoothing bandwidth that was used.
}
\details{
  The pair correlation function \eqn{g(r)} 
  is a summary of the dependence between points in a spatial point
  process. The best intuitive interpretation is the following: the probability
  \eqn{p(r)} of finding two points at locations \eqn{x} and \eqn{y}
  separated by a distance \eqn{r} is equal to
  \deqn{
    p(r) = \lambda^2 g(r) \,{\rm d}x \, {\rm d}y
  }{
    p(r) = lambda^2 * g(r) dx dy
  }
  where \eqn{\lambda}{lambda} is the intensity of the point process.
  For a completely random (uniform Poisson) process,
  \eqn{p(r) = \lambda^2 \,{\rm d}x \, {\rm d}y}{p(r) = lambda^2 dx dy}
  so \eqn{g(r) = 1}.
  Formally, the pair correlation function of a stationary point process
  is defined by 
  \deqn{
    g(r) = \frac{K'(r)}{2\pi r}
  }{
    g(r) = K'(r)/ ( 2 * pi * r) 
  }
  where \eqn{K'(r)} is the derivative of \eqn{K(r)}, the
  reduced second moment function (aka ``Ripley's \eqn{K} function'')
  of the point process. See \code{\link{Kest}} for information
  about \eqn{K(r)}.

  For a stationary Poisson process, the
  pair correlation function is identically equal to 1. Values
  \eqn{g(r) < 1} suggest inhibition between points;
  values greater than 1 suggest clustering.

  This routine computes an estimate of \eqn{g(r)}
  by kernel smoothing. 

  \itemize{
    \item
    If \code{divisor="r"} (the default), then the standard
    kernel estimator (Stoyan and Stoyan, 1994, pages 284--285)
    is used. By default, the recommendations of Stoyan and Stoyan (1994)
    are followed exactly. 
    \item
    If \code{divisor="d"} then a modified estimator is used:
    the contribution from
    an interpoint distance \eqn{d_{ij}}{d[ij]} to the
    estimate of \eqn{g(r)} is divided by \eqn{d_{ij}}{d[ij]}
    instead of dividing by \eqn{r}. This usually improves the
    bias of the estimator when \eqn{r} is close to zero.
  }

  There is also a choice of spatial edge corrections
  (which are needed to avoid bias due to edge effects
  associated with the boundary of the spatial window):

  \itemize{
    \item
    If \code{correction="translate"} or \code{correction="translation"}
    then the translation correction
    is used. For \code{divisor="r"} the translation-corrected estimate
    is given in equation (15.15), page 284 of Stoyan and Stoyan (1994).
    \item
    If \code{correction="Ripley"} then Ripley's isotropic edge correction
    is used. For \code{divisor="r"} the isotropic-corrected estimate
    is given in equation (15.18), page 285 of Stoyan and Stoyan (1994). 
    \item
    If \code{correction=c("translate", "Ripley")} then both estimates
    will be computed.
  }
  Alternatively \code{correction="all"} selects all options.
  
  The choice of smoothing kernel is controlled by the 
  argument \code{kernel} which is passed to \code{\link{density}}.
  The default is the Epanechnikov kernel, recommended by
  Stoyan and Stoyan (1994, page 285).

  The bandwidth of the smoothing kernel can be controlled by the
  argument \code{bw}. Its precise interpretation
  is explained in the documentation for \code{\link{density.default}}.
  For the Epanechnikov kernel, the argument \code{bw} is
  equivalent to \eqn{h/\sqrt{5}}{h/sqrt(5)}.

  Stoyan and Stoyan (1994, page 285) recommend using the Epanechnikov
  kernel with support \eqn{[-h,h]} chosen by the rule of thumn
  \eqn{h = c/\sqrt{\lambda}}{h = c/sqrt(lambda)},
  where \eqn{\lambda}{lambda} is the (estimated) intensity of the
  point process, and \eqn{c} is a constant in the range from 0.1 to 0.2.
  See equation (15.16).
  If \code{bw} is missing, then this rule of thumb will be applied.
  The argument \code{stoyan} determines the value of \eqn{c}.
  The smoothing bandwidth that was used in the calculation is returned
  as an attribute of the final result.
  
  The argument \code{r} is the vector of values for the
  distance \eqn{r} at which \eqn{g(r)} should be evaluated.
  There is a sensible default.
  If it is specified, \code{r} must be a vector of increasing numbers
  starting from \code{r[1] = 0}, 
  and \code{max(r)} must not exceed half the diameter of 
  the window.

  If the argument \code{domain} is given, estimation will be restricted
  to this region. That is, the estimate of 
  \eqn{g(r)} will be based on pairs of points in which the first point lies
  inside \code{domain} and the second point is unrestricted.
  The argument \code{domain}
  should be a window (object of class \code{"owin"}) or something acceptable to
  \code{\link{as.owin}}. It must be a subset of the
  window of the point pattern \code{X}.

  To compute a confidence band for the true value of the
  pair correlation function, use \code{\link{lohboot}}.

  If \code{var.approx = TRUE}, the variance of the
  estimate of the pair correlation will also be calculated using
  an analytic approximation (Illian et al, 2008, page 234)
  which is valid for stationary point processes which are not
  too clustered. This calculation is not yet implemented when
  the argument \code{domain} is given.
}
\references{
  Illian, J., Penttinen, A., Stoyan, H. and Stoyan, D. (2008)
  \emph{Statistical Analysis and Modelling of Spatial Point Patterns.}
  Wiley.

  Stoyan, D. and Stoyan, H. (1994)
  \emph{Fractals, random shapes and point fields:
  methods of geometrical statistics.}
  John Wiley and Sons.
}
\seealso{
  \code{\link{Kest}},
  \code{\link{pcf}},
  \code{\link{density}},
  \code{\link{lohboot}}.
}
\examples{
  X <- simdat
  \testonly{
    X <- X[seq(1,npoints(X), by=4)]
  }
  p <- pcf(X)
  plot(p, main="pair correlation function for X")
  # indicates inhibition at distances r < 0.3

  pd <- pcf(X, divisor="d")

  # compare estimates
  plot(p, cbind(iso, theo) ~ r, col=c("blue", "red"),
         ylim.covers=0, main="", lwd=c(2,1), lty=c(1,3), legend=FALSE)
  plot(pd, iso ~ r, col="green", lwd=2, add=TRUE)
  legend("center", col=c("blue", "green"), lty=1, lwd=2,
         legend=c("divisor=r","divisor=d"))

  # calculate approximate variance and show POINTWISE confidence bands
  pv <- pcf(X, var.approx=TRUE)
  plot(pv, cbind(iso, iso+2*sqrt(v), iso-2*sqrt(v)) ~ r)
}
\author{
  \spatstatAuthors
  and Martin Hazelton.
}
\keyword{spatial}
\keyword{nonparametric}
