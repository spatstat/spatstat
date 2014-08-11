\name{pcf.ppp}
\alias{pcf.ppp}
\title{Pair Correlation Function of Point Pattern}
\description{
  Estimates the pair correlation function of
  a point pattern using kernel methods.
}
\usage{
  \method{pcf}{ppp}(X, \dots, r = NULL, kernel="epanechnikov", bw=NULL, stoyan=0.15,
                    correction=c("translate", "Ripley"))
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
  as required.
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
  \eqn{p(r) = \lambda^2}{p(r) = lambda^2} so \eqn{g(r) = 1}.

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
  by the kernel smoothing method (Stoyan and Stoyan (1994),
  pages 284--285). By default, their recommendations
  are followed exactly.

  If \code{correction="translate"} or \code{correction="translation"}
  then the translation correction
  is used. The estimate is given in equation (15.15), page 284 of
  Stoyan and Stoyan (1994).

  If \code{correction="Ripley"} then Ripley's isotropic edge correction
  is used; the estimate is given in equation (15.18), page 285
  of Stoyan and Stoyan (1994).

  If \code{correction=c("translate", "Ripley")} then both estimates
  will be computed.

  The choice of smoothing kernel is controlled by the 
  argument \code{kernel} which is passed to \code{\link{density}}.
  The default is the Epanechnikov kernel, recommended by
  Stoyan and Stoyan (1994, page 285).

  The bandwidth of the smoothing kernel can be controlled by the
  argument \code{bw}. Its precise interpretation
  is explained in the documentation for \code{\link{density}}.
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

  The argument \code{r} is the vector of values for the
  distance \eqn{r} at which \eqn{g(r)} should be evaluated.
  There is a sensible default.
  If it is specified, \code{r} must be a vector of increasing numbers
  starting from \code{r[1] = 0}, 
  and \code{max(r)} must not exceed half the diameter of 
  the window.

  To compute a confidence band for the true value of the
  pair correlation function, use \code{\link{lohboot}}.
}
\references{
  Stoyan, D. and Stoyan, H. (1994)
  Fractals, random shapes and point fields:
  methods of geometrical statistics.
  John Wiley and Sons.
}
\seealso{
  \code{\link{Kest}},
  \code{\link{pcf}},
  \code{\link{density}},
  \code{\link{lohboot}}.
}
\examples{
  data(simdat)
  \testonly{
    simdat <- simdat[seq(1,simdat$n, by=4)]
  }
  p <- pcf(simdat)
  plot(p, main="pair correlation function for simdat")
  # indicates inhibition at distances r < 0.3
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{nonparametric}
