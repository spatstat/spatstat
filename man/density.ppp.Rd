\name{density.ppp}
\alias{density.ppp}
\title{Kernel Smoothed Intensity of Point Pattern}
\description{
  Compute a kernel smoothed intensity function from a point pattern.
}
\usage{
  \method{density}{ppp}(x, sigma, \dots,
        weights, edge=TRUE, varcov=NULL,
        at="pixels", leaveoneout=TRUE,
        adjust=1, diggle=FALSE)
}
\arguments{
  \item{x}{
    Point pattern (object of class \code{"ppp"}).
  }
  \item{sigma}{
    Standard deviation of isotropic Gaussian smoothing kernel.
    Either a numerical value, or a function that computes an
    appropriate value of \code{sigma}.
  }
  \item{weights}{
    Optional vector or matrix of weights to be attached to the points.
    May include negative values. 
  }
  \item{\dots}{
    Arguments passed to \code{\link{as.mask}} to determine
    the pixel resolution.
  }
  \item{edge}{
    Logical flag: if \code{TRUE}, apply edge correction.
  }
  \item{varcov}{
    Variance-covariance matrix of anisotropic Gaussian kernel.
    Incompatible with \code{sigma}.
  }
  \item{at}{
    String specifying whether to compute the intensity values
    at a grid of pixel locations (\code{at="pixels"}) or
    only at the points of \code{x} (\code{at="points"}).
  }
  \item{leaveoneout}{
    Logical value indicating whether to compute a leave-one-out
    estimator. Applicable only when \code{at="points"}.
  }
  \item{adjust}{
    Optional. Adjustment factor for the smoothing parameter.
  }
  \item{diggle}{
    Logical. If \code{TRUE}, use Diggle's edge correction,
    which is more accurate but slower to compute than the
    correction described under Details.
  }
}
\value{
  By default, the result is
  a pixel image (object of class \code{"im"}). 
  Pixel values are estimated intensity values,
  expressed in \dQuote{points per unit area}.

  If \code{at="points"}, the result is a numeric vector
  of length equal to the number of points in \code{x}.
  Values are estimated intensity values at the points of \code{x}.

  In either case, the return value has attributes
  \code{"sigma"} and \code{"varcov"} which report the smoothing
  bandwidth that was used.

  If \code{weights} is a matrix with more than one column, then the
  result is a list of images (if \code{at="pixels"}) or a matrix of
  numerical values (if \code{at="points"}).
  
}
\details{
  This is a method for the generic function \code{density}.

  It computes a fixed-bandwidth kernel estimate 
  (Diggle, 1985) of the intensity function of the point process
  that generated the point pattern \code{x}.
  
  By default it computes the convolution of the
  isotropic Gaussian kernel of standard deviation \code{sigma}
  with point masses at each of the data points in \code{x}.
  Anisotropic Gaussian kernels are also supported.
  Each point has unit weight, unless the argument \code{weights} is
  given (it should be a numeric vector or matrix;
  weights can be negative or zero).

  If \code{edge=TRUE}, the intensity estimate is corrected for
  edge effect bias in one of two ways:
  \itemize{
    \item If \code{diggle=FALSE} (the default) the intensity estimate is
    correted by dividing it by the convolution of the
    Gaussian kernel with the window of observation.
    Thus the intensity value at a point \eqn{u} is
    \deqn{
      \hat\lambda(u) = e(u) \sum_i k(x_i - u) w_i
    }{
      lambda(u) = e(u) sum[i] k(x[i] - u) w[i]
    }
    where \eqn{k} is the Gaussian smoothing kernel,
    \eqn{e(u)} is an edge correction factor, 
    and \eqn{w_i}{w[i]} are the weights.
    \item
    If \code{diggle=TRUE} then the method of Diggle (1985)
    is followed exactly.
    The intensity value at a point \eqn{u} is 
    \deqn{
      \hat\lambda(u) = \sum_i k(x_i - u) w_i e(x_i)
    }{
      lambda(u) = sum[i] k(x[i] - u) w[i] e(x[i])
    }
    where again \eqn{k} is the Gaussian smoothing kernel,
    \eqn{e(x_i)}{e(x[i])} is an edge correction factor, 
    and \eqn{w_i}{w[i]} are the weights.
    This computation is slightly slower but more accurate.
  }
  In both cases, the edge correction term \eqn{e(u)} is the reciprocal of the
  kernel mass inside the window:
  \deqn{
    \frac{1}{e(u)} = \int_W k(v-u) \, {\rm d}v
  }{
    1/e(u) = integral[v in W] k(v-u) dv
  }
  where \eqn{W} is the observation window.

  The smoothing kernel is determined by the arguments
  \code{sigma}, \code{varcov} and \code{adjust}.
  \itemize{
    \item if \code{sigma} is a single numerical value,
    this is taken as the standard deviation of the isotropic Gaussian
    kernel.
    \item alternatively \code{sigma} may be a function that computes
    an appropriate bandwidth for the isotropic Gaussian kernel
    from the data point pattern by calling \code{sigma(x)}.
    To perform automatic bandwidth selection using cross-validation,
    it is recommended to use the function \code{\link{bw.diggle}}.
    \item
    The smoothing kernel may be chosen to be any Gaussian
    kernel, by giving the variance-covariance matrix \code{varcov}.
    The arguments \code{sigma} and \code{varcov} are incompatible.
    \item
    Alternatively \code{sigma} may be a vector of length 2 giving the
    standard deviations of two independent Gaussian coordinates,
    thus equivalent to \code{varcov = diag(rep(sigma^2, 2))}.
    \item if neither \code{sigma} nor \code{varcov} is specified,
    an isotropic Gaussian kernel will be used, 
    with a default value of \code{sigma}
    calculated by a simple rule of thumb
    that depends only on the size of the window.
    \item
    The argument \code{adjust} makes it easy for the user to change the
    bandwidth specified by any of the rules above.
    The value of \code{sigma} will be multiplied by
    the factor \code{adjust}. The matrix \code{varcov} will be
    multiplied by \code{adjust^2}. To double the smoothing bandwidth, set
    \code{adjust=2}.
  }
  
  By default the intensity values are
  computed at every location \eqn{u} in a fine grid,
  and are returned as a pixel image.
  Computation is performed using the Fast Fourier Transform.
  Accuracy depends on the pixel resolution, controlled by the arguments
  \code{\dots} passed to \code{\link{as.mask}}.

  If \code{at="points"}, the intensity values are computed 
  to high accuracy at the points of \code{x} only. Computation is
  performed by directly evaluating and summing the Gaussian kernel
  contributions without discretising the data. The result is a numeric
  vector giving the density values.
  The intensity value at a point \eqn{x_i}{x[i]} is (if \code{diggle=FALSE})
  \deqn{
    \hat\lambda(x_i) = e(x_i) \sum_j k(x_j - x_i) w_j
  }{
    lambda(x[i]) = e(x[i]) sum[j] k(x[j] - x[i]) w[j]
  }
  or (if \code{diggle=TRUE})
  \deqn{
    \hat\lambda(x_i) = \sum_j k(x_j - x_i) w_j e(x_j)
  }{
    lambda(x[i]) = sum[j] k(x[j] - x[i]) w[j] e(x[j])
  }
  If \code{leaveoneout=TRUE} (the default), then the sum in the equation
  is taken over all \eqn{j} not equal to \eqn{i},
  so that the intensity value at a
  data point is the sum of kernel contributions from
  all \emph{other} data points.
  If \code{leaveoneout=FALSE} then the sum is taken over all \eqn{j},
  so that the intensity value at a data point includes a contribution
  from the same point.
  
  If \code{weights} is a matrix with more than one column, then the
  calculation is effectively repeated for each column of weights. The
  result is a list of images (if \code{at="pixels"}) or a matrix of
  numerical values (if \code{at="points"}).
  
  To select the bandwidth \code{sigma} automatically by
  cross-validation, use \code{\link{bw.diggle}}.
  
  To perform spatial interpolation of values that were observed
  at the points of a point pattern, use \code{\link{smooth.ppp}}.

  For adaptive nonparametric estimation, see
  \code{\link{adaptive.density}}.
  For data sharpening, see \code{\link{sharpen.ppp}}.

  To compute a relative risk surface or probability map for
  two (or more) types of points, use \code{\link{relrisk}}.
}
\seealso{
  \code{\link{bw.diggle}},
  \code{\link{smooth.ppp}},
  \code{\link{sharpen.ppp}},
  \code{\link{adaptive.density}},
  \code{\link{relrisk}},
  \code{\link{ppp.object}},
  \code{\link{im.object}}
}
\note{
  This function is often misunderstood.

  The result of \code{density.ppp} is not a spatial smoothing 
  of the marks or weights attached to the point pattern.
  To perform spatial interpolation of values that were observed
  at the points of a point pattern, use \code{\link{smooth.ppp}}.

  The result of \code{density.ppp} is not a probability density.
  It is an estimate of the \emph{intensity function} of the
  point process that generated the point pattern data.
  Intensity is the expected number of random points
  per unit area.
  The units of intensity are \dQuote{points per unit area}.
  Intensity is usually a function of spatial location,
  and it is this function which is estimated by \code{density.ppp}.
  The integral of the intensity function over a spatial region gives the
  expected number of points falling in this region.

  Inspecting an estimate of the intensity function is usually the
  first step in exploring a spatial point pattern dataset.
  For more explanation, see the workshop notes (Baddeley, 2008)
  or Diggle (2003).

  If you have two (or more) types of points, and you want a
  probability map or relative risk surface (the spatially-varying
  probability of a given type), use \code{\link{relrisk}}.
}
\examples{
  data(cells)
  if(interactive()) {
    opa <- par(mfrow=c(1,2))
    plot(density(cells, 0.05))
    plot(density(cells, 0.05, diggle=TRUE))
    par(opa)
    v <- diag(c(0.05, 0.07)^2)
    plot(density(cells, varcov=v))
  }
  \testonly{
     Z <- density(cells, 0.05)
     Z <- density(cells, 0.05, diggle=TRUE)
     Z <- density(cells, varcov=diag(c(0.05^2, 0.07^2)))
     Z <- density(cells, 0.05, weights=data.frame(a=1:42,b=42:1))
  }
  # automatic bandwidth selection
  plot(density(cells, sigma=bw.diggle(cells)))
  # equivalent:
  plot(density(cells, bw.diggle))
  # evaluate intensity at points
  density(cells, 0.05, at="points")
}
\references{
  Baddeley, A. (2010) Analysing spatial point patterns in R.
  Workshop notes. CSIRO online technical publication.
  URL: \code{www.csiro.au/resources/pf16h.html}

  Diggle, P.J. (1985)
  A kernel method for smoothing point process data.
  \emph{Applied Statistics} (Journal of the Royal Statistical Society,
  Series C) \bold{34} (1985) 138--147.

  Diggle, P.J. (2003)
  \emph{Statistical analysis of spatial point patterns},
  Second edition. Arnold.
}

\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{methods}
\keyword{smooth}
