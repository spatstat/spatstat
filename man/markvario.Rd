\name{markvario}
\alias{markvario}
\title{Mark Variogram}
\description{
  Estimate the mark variogram of a marked point pattern.
}
\usage{
markvario(X, correction = c("isotropic", "Ripley", "translate"),
r = NULL, method = "density", ..., normalise=FALSE)
}
\arguments{
  \item{X}{The observed point pattern.
    An object of class \code{"ppp"} or something acceptable to
    \code{\link{as.ppp}}. It must have marks which are numeric.
  }
  \item{correction}{
    A character vector containing any selection of the
    options \code{"isotropic"}, \code{"Ripley"} or \code{"translate"}.
    It specifies the edge correction(s) to be applied.
  }
  \item{r}{numeric vector. The values of the argument \eqn{r}
    at which the mark variogram
    \eqn{\gamma(r)}{gamma(r)}
    should be evaluated.
    There is a sensible default.
  }
  \item{method}{
    A character vector indicating the user's choice of
    density estimation technique to be used. Options are
    \code{"density"}, 
    \code{"loess"},
    \code{"sm"} and \code{"smrep"}.
  }
  \item{\dots}{
    Other arguments passed to \code{\link{markcorr}},
    or passed to the density estimation routine
    (\code{\link{density}}, \code{\link{loess}} or \code{sm.density})
    selected by \code{method}.
  }
  \item{normalise}{If \code{TRUE}, normalise the variogram by
    dividing it by the estimated mark variance.
  }
}
\details{
  The mark variogram \eqn{\gamma(r)}{gamma(r)}
  of a marked point process \eqn{X}
  is a measure of the dependence between the marks of two 
  points of the process a distance \eqn{r} apart.
  It is informally defined as
  \deqn{
    \gamma(r) = E[\frac 1 2 (M_1 - M_2)^2]
  }{
    gamma(r) = E[(1/2) * (M1 - M2)^2 ]
  }
  where \eqn{E[ ]} denotes expectation and \eqn{M_1,M_2}{M1,M2}
  are the marks attached to two points of the process
  a distance \eqn{r} apart.

  The mark variogram of a marked point process is analogous,
  but \bold{not equivalent}, to the variogram of a random field
  in geostatistics. See Waelder and Stoyan (1996).
}
\value{
  An object of class \code{"fv"} (see \code{\link{fv.object}}).
  
  Essentially a data frame containing numeric columns 
  \item{r}{the values of the argument \eqn{r} 
    at which the mark variogram \eqn{\gamma(r)}{gamma(r)}
    has been  estimated
  }
  \item{theo}{the theoretical value of \eqn{\gamma(r)}{gamma(r)}
    when the marks attached to different points are independent;
    equal to the sample variance of the marks
  }
  together with a column or columns named 
  \code{"iso"} and/or \code{"trans"},
  according to the selected edge corrections. These columns contain
  estimates of the function \eqn{\gamma(r)}{gamma(r)}
  obtained by the edge corrections named.
}
\references{
  Cressie, N.A.C. (1991)
  \emph{Statistics for spatial data}.
  John Wiley and Sons, 1991.
  
  Mase, S. (1996) 
  The threshold method for estimating annual rainfall.
  \emph{Annals of the Institute of Statistical Mathematics}
  \bold{48} (1996) 201-213.

  Waelder, O. and Stoyan, D. (1996)
  On variograms in point process statistics.
  \emph{Biometrical Journal} \bold{38} (1996) 895-905.
}
\seealso{
  Mark correlation function \code{\link{markcorr}} for numeric marks.

  Mark connection function \code{\link{markconnect}} and 
  multitype K-functions \code{\link{Kcross}}, \code{\link{Kdot}}
  for factor-valued marks.
}
\examples{
    # Longleaf Pine data
    # marks represent tree diameter
    data(longleaf)
    # Subset of this large pattern
    swcorner <- owin(c(0,100),c(0,100))
    sub <- longleaf[ , swcorner]
    # mark correlation function
    mv <- markvario(sub)
    plot(mv)
}
\author{\adrian
  
  
  and \rolf
  
}
\keyword{spatial}
\keyword{nonparametric}
