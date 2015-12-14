\name{pcfdot}
\alias{pcfdot}
\title{Multitype pair correlation function (i-to-any)}
\description{
  Calculates an estimate of the multitype pair correlation function
  (from points of type \code{i} to points of any type)
  for a multitype point pattern.
}
\usage{
  pcfdot(X, i, ..., r = NULL,
         kernel = "epanechnikov", bw = NULL, stoyan = 0.15,
         correction = c("isotropic", "Ripley", "translate"),
         divisor = c("r", "d"))
}
\arguments{
  \item{X}{The observed point pattern, 
    from which an estimate of the dot-type pair correlation function
    \eqn{g_{i\bullet}(r)}{gdot[i](r)} will be computed.
    It must be a multitype point pattern (a marked point pattern
    whose marks are a factor). 
  }
  \item{i}{The type (mark value)
    of the points in \code{X} from which distances are measured.
    A character string (or something that will be converted to a
    character string).
    Defaults to the first level of \code{marks(X)}.
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
    Coefficient for default bandwidth rule; see Details.
  }
  \item{correction}{
    Choice of edge correction.
  }
  \item{divisor}{
    Choice of divisor in the estimation formula:
    either \code{"r"} (the default) or \code{"d"}. See Details.
  }
}
\details{
  This is a generalisation of the pair correlation function \code{\link{pcf}}
  to multitype point patterns.

  For two locations \eqn{x} and \eqn{y} separated by a nonzero
  distance \eqn{r},
  the probability \eqn{p(r)} of finding a point of type \eqn{i} at location
  \eqn{x} and a point of any type at location \eqn{y} is 
  \deqn{
    p(r) = \lambda_i \lambda g_{i\bullet}(r) \,{\rm d}x \, {\rm d}y
  }{
    p(r) = lambda[i] * lambda * gdot[i](r) dx dy
  }
  where \eqn{\lambda}{lambda} is the intensity of all points,
  and \eqn{\lambda_i}{lambda[i]} is the intensity of the points
  of type \eqn{i}. 
  For a completely random Poisson marked point process,
  \eqn{p(r) = \lambda_i \lambda}{p(r) = lambda[i] * lambda}
  so \eqn{g_{i\bullet}(r) = 1}{gdot[i](r) = 1}.
  
  For a stationary multitype point process, the
  type-\code{i}-to-any-type pair correlation
  function between marks \eqn{i} and \eqn{j} is formally defined as
  \deqn{
    g_{i\bullet}(r) = \frac{K_{i\bullet}^\prime(r)}{2\pi r}
  }{
    g(r) = Kdot[i]'(r)/ ( 2 * pi * r) 
  }
  where \eqn{K_{i\bullet}^\prime}{Kdot[i]'(r)} is the derivative of
  the type-\code{i}-to-any-type \eqn{K} function
  \eqn{K_{i\bullet}(r)}{Kdot[i](r)}.
  of the point process. See \code{\link{Kdot}} for information
  about   \eqn{K_{i\bullet}(r)}{Kdot[i](r)}.

  The command \code{pcfdot} computes a kernel estimate of
  the multitype pair correlation function from points of type \eqn{i}
  to points of any type.

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
  285). That is,
  \eqn{h = c/\sqrt{\lambda}}{h = c/sqrt(lambda)},
  where \eqn{\lambda}{lambda} is the (estimated) intensity of the
  unmarked point process, 
  and \eqn{c} is a constant in the range from 0.1 to 0.2.
  The argument \code{stoyan} determines the value of \eqn{c}.

  The companion function \code{\link{pcfcross}} computes the
  corresponding analogue of \code{\link{Kcross}}.
}

\value{
  An object of class \code{"fv"}, see \code{\link{fv.object}},
  which can be plotted directly using \code{\link{plot.fv}}.

  Essentially a data frame containing columns
  \item{r}{the vector of values of the argument \eqn{r} 
    at which the function \eqn{g_{i\bullet}}{gdot[i]} has been  estimated
  }
  \item{theo}{the theoretical value \eqn{g_{i\bullet}(r) = 1}{gdot[i](r) = r}
    for independent marks.
  }
  together with columns named 
  \code{"border"}, \code{"bord.modif"},
  \code{"iso"} and/or \code{"trans"},
  according to the selected edge corrections. These columns contain
  estimates of the function \eqn{g_{i,j}}{g[i,j]}
  obtained by the edge corrections named.
}
\seealso{
  Mark connection function \code{\link{markconnect}}.

  Multitype pair correlation \code{\link{pcfcross}}, \code{\link{pcfmulti}}.
  
  Pair correlation \code{\link{pcf}},\code{\link{pcf.ppp}}.
  
  \code{\link{Kdot}}
}
\examples{
 data(amacrine)
 p <- pcfdot(amacrine, "on")
 p <- pcfdot(amacrine, "on", stoyan=0.1)
 plot(p)
}
\author{Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  and Rolf Turner \email{r.turner@auckland.ac.nz}
  
}
\keyword{spatial}
\keyword{nonparametric}
