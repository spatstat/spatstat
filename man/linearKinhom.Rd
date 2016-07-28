\name{linearKinhom}
\alias{linearKinhom}
\title{
  Inhomogeneous Linear K Function
}
\description{
  Computes an estimate of the inhomogeneous linear \eqn{K} function
  for a point pattern on a linear network.
}
\usage{
linearKinhom(X, lambda=NULL, r=NULL, ..., correction="Ang",
             normalise=TRUE, normpower=1)
}
\arguments{
  \item{X}{
    Point pattern on linear network (object of class \code{"lpp"}).
  }
  \item{lambda}{
    Intensity values for the point pattern. Either a numeric vector,
    a \code{function}, a pixel image
    (object of class \code{"im"} or \code{"linim"}) or
    a fitted point process model (object of class \code{"ppm"}
    or \code{"lppm"}).
  }
  \item{r}{
    Optional. Numeric vector of values of the function argument \eqn{r}.
    There is a sensible default.
  }
  \item{\dots}{
    Ignored.
  }
  \item{correction}{
    Geometry correction.
    Either \code{"none"} or \code{"Ang"}. See Details.
  }
  \item{normalise}{
    Logical. If \code{TRUE} (the default), the denominator of the estimator is 
    data-dependent (equal to the sum of the reciprocal intensities at the data
    points), which reduces the sampling variability.
    If \code{FALSE}, the denominator is the length of the network.
  }
  \item{normpower}{
    Integer (usually either 1 or 2).
    Normalisation power. See Details.
  }
}
\details{
  This command computes the inhomogeneous version of the 
  linear \eqn{K} function from point pattern data on a linear network.

  If \code{lambda = NULL} the result is equivalent to the
  homogeneous \eqn{K} function \code{\link{linearK}}.
  If \code{lambda} is given, then it is expected to provide estimated values
  of the intensity of the point process at each point of \code{X}. 
  The argument \code{lambda} may be a numeric vector (of length equal to
  the number of points in \code{X}), or a \code{function(x,y)} that will be
  evaluated at the points of \code{X} to yield numeric values, 
  or a pixel image (object of class \code{"im"}) or a fitted point 
  process model (object of class \code{"ppm"} or \code{"lppm"}).

  If \code{lambda} is a fitted point process model,
  the default behaviour is to update the model by re-fitting it to
  the data, before computing the fitted intensity.
  This can be disabled by setting \code{update=FALSE}.
  
  If \code{correction="none"}, the calculations do not include
  any correction for the geometry of the linear network.
  If \code{correction="Ang"}, the pair counts are weighted using
  Ang's correction (Ang, 2010).

  Each estimate is initially computed as 
  \deqn{
    \widehat K_{\rm inhom}(r) = \frac{1}{\mbox{length}(L)}
    \sum_i \sum_j \frac{1\{d_{ij} \le r\}
      e(x_i,x_j)}{\lambda(x_i)\lambda(x_j)}
  }{
    K^inhom(r)= (1/length(L)) sum[i] sum[j] 1(d[i,j] <= r) * 
    e(x[i],x[j])/(lambda(x[i]) * lambda(x[j]))
  }
  where \code{L} is the linear network,
  \eqn{d_{ij}}{d[i,j]} is the distance between points
  \eqn{x_i}{x[i]} and \eqn{x_j}{x[j]}, and
  \eqn{e(x_i,x_j)}{e(x[i],x[j])} is a weight.
  If \code{correction="none"} then this weight is equal to 1,
  while if  \code{correction="Ang"} the weight is
  \eqn{e(x_i,x_j,r) = 1/m(x_i, d_{ij})}{e(x[i],x[j],r) = 1/m(x[i],d[i,j])}
  where \eqn{m(u,t)} is the number of locations on the network that lie
  exactly \eqn{t} units distant from location \eqn{u} by the shortest
  path.

  If \code{normalise=TRUE} (the default), then the estimates
  described above
  are multiplied by \eqn{c^{\mbox{normpower}}}{c^normpower} where 
  \eqn{
    c = \mbox{length}(L)/\sum (1/\lambda(x_i)).
  }{
    c = length(L)/sum[i] (1/lambda(x[i])).
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
}
\value{
  Function value table (object of class \code{"fv"}).
}
\author{
  Ang Qi Wei \email{aqw07398@hotmail.com} and
  \adrian
}
\references{
  Ang, Q.W. (2010) Statistical methodology for spatial point patterns
  on a linear network. MSc thesis, University of Western Australia.
  
  Ang, Q.W., Baddeley, A. and Nair, G. (2012)
  Geometrically corrected second-order analysis of 
  events on a linear network, with applications to
  ecology and criminology.
  \emph{Scandinavian Journal of Statistics} \bold{39}, 591--617.
}
\seealso{
  \code{\link{lpp}}
}
\examples{
  data(simplenet)
  X <- rpoislpp(5, simplenet)
  fit <- lppm(X, ~x)
  K <- linearKinhom(X, lambda=fit)
  plot(K)
}
\keyword{spatial}
\keyword{nonparametric}


