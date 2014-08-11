\name{pcfmulti}
\alias{pcfmulti}
\title{
  Marked pair correlation function
}
\description{
For a marked point pattern, 
estimate the multitype pair correlation function
using kernel methods.
}
\usage{
   pcfmulti(X, I, J, ..., r = NULL,
            kernel = "epanechnikov", bw = NULL, stoyan = 0.15,
            correction = c("translate", "Ripley"),
            divisor = c("r", "d"),
            Iname = "points satisfying condition I",
            Jname = "points satisfying condition J")
}
\arguments{
  \item{X}{The observed point pattern, 
    from which an estimate of the cross-type pair correlation function
    \eqn{g_{ij}(r)}{g[i,j](r)} will be computed.
    It must be a multitype point pattern (a marked point pattern
    whose marks are a factor). 
  }
  \item{I}{Subset index specifying the points of \code{X}
    from which distances are measured. 
  }
  \item{J}{Subset index specifying the points in \code{X} to which
    distances are measured. 
  }
  \item{\dots}{
    Ignored.
  }
  \item{r}{
    Vector of values for the argument \eqn{r} at which \eqn{g(r)} 
    should be evaluated. There is a sensible default.
  }
  \item{kernel}{
    Choice of smoothing kernel,
    passed to \code{\link{density.default}}.
  }
  \item{bw}{
    Bandwidth for smoothing kernel, 
    passed to \code{\link{density.default}}.
  }
  \item{stoyan}{
    Coefficient for default bandwidth rule.
  }
  \item{correction}{
    Choice of edge correction.
  }
  \item{divisor}{
    Choice of divisor in the estimation formula:
    either \code{"r"} (the default) or \code{"d"}. 
  }
  \item{Iname,Jname}{
    Optional. Character strings describing the members of
    the subsets \code{I} and \code{J}.
  }
}
\details{
  This is a generalisation of \code{\link{pcfcross}}
  to arbitrary collections of points.

  The algorithm measures the distance from each data point
  in subset \code{I} to each data point in subset \code{J},
  excluding identical pairs of points. The distances are
  kernel-smoothed and renormalised to form a pair correlation
  function.

  \itemize{
    \item
    If \code{divisor="r"} (the default), then the multitype
    counterpart of the standard
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
  \code{correction="translate"} is the Ohser-Stoyan translation
  correction, and \code{correction="isotropic"} or \code{"Ripley"}
  is Ripley's isotropic correction.  

  The arguments \code{I} and \code{J} specify two subsets of the
  point pattern \code{X}. They may be any type of subset indices, for example,
  logical vectors of length equal to \code{npoints(X)},
  or integer vectors with entries in the range 1 to
  \code{npoints(X)}, or negative integer vectors.

  Alternatively, \code{I} and \code{J} may be \bold{functions}
  that will be applied to the point pattern \code{X} to obtain
  index vectors. If \code{I} is a function, then evaluating
  \code{I(X)} should yield a valid subset index. This option
  is useful when generating simulation envelopes using
  \code{\link{envelope}}.

  The choice of smoothing kernel is controlled by the 
  argument \code{kernel} which is passed to \code{\link{density}}.
  The default is the Epanechnikov kernel.

  The bandwidth of the smoothing kernel can be controlled by the
  argument \code{bw}. Its precise interpretation
  is explained in the documentation for \code{\link{density.default}}.
  For the Epanechnikov kernel with support \eqn{[-h,h]},
  the argument \code{bw} is equivalent to \eqn{h/\sqrt{5}}{h/sqrt(5)}.

  If \code{bw} is not specified, the default bandwidth
  is determined by Stoyan's rule of thumb (Stoyan and Stoyan, 1994, page
  285) applied to the points of type \code{j}. That is,
  \eqn{h = c/\sqrt{\lambda}}{h = c/sqrt(lambda)},
  where \eqn{\lambda}{lambda} is the (estimated) intensity of the
  point process of type \code{j},
  and \eqn{c} is a constant in the range from 0.1 to 0.2.
  The argument \code{stoyan} determines the value of \eqn{c}.
}
\value{
  An object of class \code{"fv"}.
}
\seealso{
  \code{\link{pcfcross}}, 
  \code{\link{pcfdot}}, 
  \code{\link{pcf.ppp}}.
}
\examples{
  adult <- (marks(longleaf) >= 30)
  juvenile <- !adult
  p <- pcfmulti(longleaf, adult, juvenile)
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{nonparametric}

