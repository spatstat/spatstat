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
bw.abram(X, h0, at=c("points", "pixels"),
         \dots, hp = NULL, pilot = NULL, trim=5)
}
\arguments{
  \item{X}{
    A point pattern (object of class \code{"ppp"})
    for which the variable bandwidths should be computed.
  }
  \item{h0}{
    A scalar value giving the global smoothing bandwidth
    in the same units as the coordinates of \code{X}.
  }
  \item{at}{
    Character string (partially matched) specifying whether
    to compute bandwidth values at the points of \code{X}
    (\code{at="points"}, the default) or to compute bandwidths
    at every pixel in a fine pixel grid (\code{at="pixels"}).
  }
  \item{\dots}{
    Additional arguments passed to
    either \code{\link[spatstat]{density.ppp}}
    and \code{\link{as.im}} to control the pixel resolution
    and the type of smoothing.
  }
  \item{hp}{
    Optional. A scalar pilot bandwidth, for fixed-bandwidth estimation
    of the pilot density. Ignored if \code{pilot} is a pixel image
    (object of class \code{"im"}); see below.
  }
  \item{pilot}{
    Optional. Specification of a pilot density.
    If \code{pilot=NULL} (the default) the pilot density
    is computed using a fixed-bandwidth estimate based
    on \code{X} and \code{hp}.
    If \code{pilot} is a point pattern, the pilot density is 
    is computed using a fixed-bandwidth estimate based
    on \code{pilot} and \code{hp}.
    If \code{pilot} is a pixel image (object of class \code{"im"}),
    the function assumes this to be the (possibly unnormalised)
    pilot density and \code{hp} is ignored.
    In both the latter two cases the function performs a basic check
    to ensure the spatial domains are the same as that of \code{X}.
  }
  \item{trim}{
    A trimming value required to curb excessively large bandwidths.
    See Details. The default is sensible in most cases.
  }
}
\details{
  This function is based on the results of Abramson (1982)
  and Hall and Marron (1988) for computation of adaptive bandwidths
  attached to each observation in \code{X}.

  Under this framework, we have \cr

  \eqn{h(u) =} \code{h0} * min[ \eqn{\tilde{f}(u)^{-1/2} / \gamma}, \code{trim} ],

  where \eqn{\tilde{f}(u)} is a pilot estimate of the spatially varying
  density. The variable bandwidths are rescaled by \eqn{\gamma}, the
  geometric mean of the \eqn{\tilde{f}(u)^{-1/2}} terms evaluated at the
  data; this allows the global bandwidth \code{h0} to be considered on
  the same scale as a corresponding fixed bandwidth. The trimming value
  \code{trim} has the same interpretation as the required `clipping' of
  the pilot density at some small nominal value (see Hall & Marron,
  1988), to necessarily prevent occassionally massive bandwidths (which
  can occur at very isolated observations).

  The pilot density may be directly supplied to the \code{pilot} argument
  as a pixel image; otherwise it is computed as a fixed-bandwidth kernel
  density estimate using the pilot bandwidth \code{hp}. This calculation
  may either use the same observations as in \code{X}, or a different
  pattern on the same spatial domain (performed when \code{pilot} is
  passed a spatial point pattern object).
}
\value{
  Either a numeric vector of length \code{npoints(X)}
  giving the Abramson bandwidth for each point
  (when \code{returnim = FALSE}, the default),
  or the entire pixel \code{\link[spatstat]{im}}age
  of the Abramson bandwidths over the relevant spatial domain
  (when \code{returnim = TRUE}).
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
}
\keyword{spatial}
\keyword{nonparametric}

