\name{bw.abram}
\alias{bw.abram}
\title{
  Abramson's Adaptive Bandwidths
}
\description{
  Computes adaptive smoothing bandwidths
  according to the inverse-square-root rule of Abramson (1982).
}
\usage{
bw.abram(X, h0, 
         \dots,
         at=c("points", "pixels"),
         hp = h0, pilot = NULL, trim=5, smoother=density.ppp)
}
\arguments{
  \item{X}{
    A point pattern (object of class \code{"ppp"})
    for which the variable bandwidths should be computed.
  }
  \item{h0}{
    A scalar value giving the global smoothing bandwidth
    in the same units as the coordinates of \code{X}.
    The default is \code{h0=\link{bw.ppl}(X)}.
  }
  \item{\dots}{
    Additional arguments passed to
    \code{\link{as.im}} to control the pixel resolution,
    or passed to \code{\link[spatstat]{density.ppp}} or \code{smoother}
    to control the type of smoothing, when computing the pilot estimate.
  }
  \item{at}{
    Character string (partially matched) specifying whether
    to compute bandwidth values at the points of \code{X}
    (\code{at="points"}, the default) or to compute bandwidths
    at every pixel in a fine pixel grid (\code{at="pixels"}).
  }
  \item{hp}{
    Optional. A scalar pilot bandwidth, used for estimation
    of the pilot density if required. Ignored if \code{pilot} is a pixel image
    (object of class \code{"im"}); see below.
  }
  \item{pilot}{
    Optional. Specification of a pilot density
    (possibly unnormalised).
    If \code{pilot=NULL} the pilot density
    is computed by applying fixed-bandwidth density estimation to
    \code{X} using bandwidth \code{hp}.
    If \code{pilot} is a point pattern, the pilot density is 
    is computed using a fixed-bandwidth estimate based
    on \code{pilot} and \code{hp}.
    If \code{pilot} is a pixel image (object of class \code{"im"}),
    this is taken to be the (possibly unnormalised)
    pilot density, and \code{hp} is ignored.
  }
  \item{trim}{
    A trimming value required to curb excessively large bandwidths.
    See Details. The default is sensible in most cases.
  }
  \item{smoother}{
    Smoother for the pilot.
    A function or character string, specifying the function
    to be used to compute the pilot estimate when
    \code{pilot} is \code{NULL} or is a point pattern.
  }
}
\details{
  This function computes adaptive smoothing bandwidths
  using the methods of Abramson (1982) and Hall and Marron (1988).

  If \code{at="points"} (the default) a smoothing bandwidth is
  computed for each point in the pattern \code{X}. Alternatively if
  \code{at="pixels"} a smoothing bandwidth is computed for
  each spatial location in a pixel grid.

  Under the Abramson-Hall-Marron rule, the bandwidth at location \eqn{u} is
  \deqn{
    h(u) = \mbox{\texttt{h0}}
    * \mbox{min}[ \frac{\tilde{f}(u)^{-1/2}}{\gamma}, \mbox{\texttt{trim}} ]
  }{
    h(u) = h0 * min(\tilde{f}(u)^{-1/2}/\gamma, trim)
  }
  where \eqn{\tilde{f}(u)} is a pilot estimate of the spatially varying
  probability density. The variable bandwidths are rescaled by \eqn{\gamma}, the
  geometric mean of the \eqn{\tilde{f}(u)^{-1/2}} terms evaluated at the
  data; this allows the global bandwidth \code{h0} to be considered on
  the same scale as a corresponding fixed bandwidth. The trimming value
  \code{trim} has the same interpretation as the required `clipping' of
  the pilot density at some small nominal value (see Hall and Marron,
  1988), to necessarily prevent extreme bandwidths (which
  can occur at very isolated observations).

  The pilot density or intensity is determined as follows:
  \itemize{
    \item If \code{pilot} is a pixel image, this is taken as the
    pilot density or intensity.
    \item If \code{pilot} is \code{NULL}, then the pilot intensity is
    computed as a fixed-bandwidth kernel
    intensity estimate using \code{\link{density.ppp}} applied to
    the data pattern \code{X} using
    the pilot bandwidth \code{hp}.
    \item
    If \code{pilot} is a different point pattern on the same spatial
    domain as \code{X}, then the pilot intensity is
    computed as a fixed-bandwidth kernel
    intensity estimate using \code{\link{density.ppp}} applied to
    \code{pilot} using the pilot bandwidth \code{hp}.
  }
  In each case the pilot density or intensity is renormalised to
  become a probability density, and then the Abramson rule is applied.

  Instead of calculating the pilot as a fixed-bandwidth density
  estimate, the user can specify another density estimation procedure
  using the argument \code{smoother}. This should be either a function
  or the character string name of a function. It will replace
  \code{\link{density.ppp}} as the function used to calculate the
  pilot estimate. The pilot estimate will be computed as
  \code{smoother(X, sigma=hp, ...)} if \code{pilot} is \code{NULL},
  or \code{smoother(pilot, sigma=hp, ...)} if \code{pilot} is a point
  pattern. If \code{smoother} does not recognise the argument name
  \code{sigma} for the smoothing bandwidth, then \code{hp} is effectively
  ignored, as shown in the Examples.  
}
\value{
  Either a numeric vector of length \code{npoints(X)}
  giving the Abramson bandwidth for each point
  (when \code{at = "points"}, the default),
  or the entire pixel \code{\link[spatstat]{im}}age
  of the Abramson bandwidths over the relevant spatial domain
  (when \code{at = "pixels"}).
}
\references{
  Abramson, I. (1982)
  On bandwidth variation in kernel estimates --- a square root law.
  \emph{Annals of Statistics}, \bold{10}(4), 1217-1223.\cr

  Davies, T.M. and Baddeley, A. (2018)
  Fast computation of spatially adaptive kernel estimates.
  \emph{Statistics and Computing}, \bold{28}(4), 937-956.\cr

  Davies, T.M., Marshall, J.C., and Hazelton, M.L. (2018)
  Tutorial on kernel estimation of continuous spatial
  and spatiotemporal relative risk.
  \emph{Statistics in Medicine}, \bold{37}(7), 1191-1221.\cr

  Hall, P. and Marron, J.S. (1988)
  Variable window width kernel density estimates of probability
  densities.
  \emph{Probability Theory and Related Fields}, \bold{80}, 37-49.\cr

  Silverman, B.W. (1986)
  \emph{Density Estimation for Statistics and Data Analysis}.
  Chapman and Hall, New York.
}
\author{
  Tilman M. Davies. Adapted by \adrian.
}
\examples{
# 'ch' just 58 laryngeal cancer cases
ch <- split(chorley)[[1]]

h <- bw.abram(ch,h0=1,hp=0.7)
length(h)
summary(h)
if(interactive()) hist(h)

# calculate pilot based on all 1036 observations
h.pool <- bw.abram(ch,h0=1,hp=0.7,pilot=chorley)
length(h.pool)
summary(h.pool)
if(interactive()) hist(h.pool)

# get full image used for 'h' above
him <- bw.abram(ch,h0=1,hp=0.7,at="pixels")
plot(him);points(ch,col="grey")

# use Voronoi-Dirichlet pilot ('hp' is ignored)
hvo <- bw.abram(ch, h0=1, smoother=densityVoronoi)
}
\keyword{spatial}
\keyword{nonparametric}

