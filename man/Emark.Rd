\name{Emark}
\alias{Emark}
\alias{Vmark}
\title{
  Diagnostics for random marking
}
\description{
  Estimate the summary functions \eqn{E(r)} and \eqn{V(r)} for 
  a marked point pattern, proposed by Schlather et al (2004) as diagnostics 
  for dependence between the points and the marks.
}
\usage{
Emark(X, r=NULL,
         correction=c("isotropic", "Ripley", "translate"),
         method="density", \dots, normalise=FALSE)
Vmark(X, r=NULL,
         correction=c("isotropic", "Ripley", "translate"),
         method="density", \dots, normalise=FALSE)
}
\arguments{
  \item{X}{The observed point pattern.
    An object of class \code{"ppp"} or something acceptable to
    \code{\link{as.ppp}}. The pattern should have numeric marks.
  }
  \item{r}{Optional. Numeric vector. The values of the argument \eqn{r}
    at which the function \eqn{E(r)} or \eqn{V(r)} should be evaluated.
    There is a sensible default.
  }
  \item{correction}{
    A character vector containing any selection of the
    options \code{"isotropic"}, \code{"Ripley"} or \code{"translate"}.
    It specifies the edge correction(s) to be applied.
  }
  \item{method}{
    A character vector indicating the user's choice of
    density estimation technique to be used. Options are
    \code{"density"}, 
    \code{"loess"},
    \code{"sm"} and \code{"smrep"}.
  }
  \item{\dots}{
    Arguments passed to the density estimation routine
    (\code{\link{density}}, \code{\link{loess}} or \code{sm.density})
    selected by \code{method}.
  }
  \item{normalise}{
    If\code{TRUE}, normalise the estimate of \eqn{E(r)} or \eqn{V(r)}
    so that it would have value equal to 1 if the marks are independent
    of the points.
  }
}
\value{
  If \code{marks(X)} is a numeric vector, the result is 
  an object of class \code{"fv"} (see \code{\link{fv.object}}).
  If \code{marks(X)} is a data frame, the result is
  a list of objects of class \code{"fv"}, one for each column of marks.

  An object of class \code{"fv"} is essentially
  a data frame containing numeric columns 
  \item{r}{the values of the argument \eqn{r} 
    at which the function \eqn{E(r)} or \eqn{V(r)}
    has been estimated
  }
  \item{theo}{the theoretical, constant value of \eqn{E(r)} or \eqn{V(r)}
    when the marks attached to different points are independent
  }
  together with a column or columns named 
  \code{"iso"} and/or \code{"trans"},
  according to the selected edge corrections. These columns contain
  estimates of the function \eqn{E(r)} or \eqn{V(r)}
  obtained by the edge corrections named.
}
\details{
  For a marked point process,
  Schlather et al (2004) defined the functions
  \eqn{E(r)} and \eqn{V(r)} to be the conditional mean
  and conditional variance of the mark attached to a
  typical random point, given that there exists another random
  point at a distance \eqn{r} away from it.

  More formally,
  \deqn{
    E(r) = E_{0u}[M(0)]
  }{
    E(r) = E[0u] M(0)
  }
  and
  \deqn{
    V(r) = E_{0u}[(M(0)-E(u))^2]
  }{
    V(r) = E[0u]((M(0)-E(u))^2)
  }
  where \eqn{E_{0u}}{E[0u]} denotes the conditional expectation
  given that there are points of the process at the locations
  \eqn{0} and \eqn{u} separated by a distance \eqn{r},
  and where \eqn{M(0)} denotes the mark attached to the point \eqn{0}. 

  These functions may serve as diagnostics for dependence
  between the points and the marks. If the points and marks are
  independent, then \eqn{E(r)} and \eqn{V(r)} should be
  constant (not depending on \eqn{r}). See Schlather et al (2004).

  The argument \code{X} must be a point pattern (object of class
  \code{"ppp"}) or any data that are acceptable to \code{\link{as.ppp}}.
  It must be a marked point pattern with numeric marks.

  The argument \code{r} is the vector of values for the
  distance \eqn{r} at which \eqn{k_f(r)}{k[f](r)} is estimated.

  This algorithm assumes that \code{X} can be treated
  as a realisation of a stationary (spatially homogeneous) 
  random spatial point process in the plane, observed through
  a bounded window.
  The window (which is specified in \code{X} as \code{Window(X)})
  may have arbitrary shape.

  Biases due to edge effects are
  treated in the same manner as in \code{\link{Kest}}.
  The edge corrections implemented here are
  \describe{
    \item{isotropic/Ripley}{Ripley's isotropic correction
      (see Ripley, 1988; Ohser, 1983).
      This is implemented only for rectangular and polygonal windows
      (not for binary masks).
    }
    \item{translate}{Translation correction (Ohser, 1983).
      Implemented for all window geometries, but slow for
      complex windows. 
    }
  }
  Note that the estimator assumes the process is stationary (spatially
  homogeneous). 

  The numerator and denominator of the mark correlation function
  (in the expression above) are estimated using density estimation
  techniques. The user can choose between
  \describe{
    \item{\code{"density"}}{
      which uses the standard kernel
      density estimation routine \code{\link{density}}, and
      works only for evenly-spaced \code{r} values;
    }
    \item{\code{"loess"}}{
      which uses the function \code{loess} in the
      package \pkg{modreg};
    }
    \item{\code{"sm"}}{
      which uses the function \code{sm.density} in the
      package \pkg{sm} and is extremely slow;
    }
    \item{\code{"smrep"}}{
      which uses the function \code{sm.density} in the
      package \pkg{sm} and is relatively fast, but may require manual
      control of the smoothing parameter \code{hmult}.
    }
  }
}
\references{
  Schlather, M. and Ribeiro, P. and Diggle, P. (2004)
  Detecting dependence between marks and locations of
  marked point processes.
  \emph{Journal of the Royal Statistical Society, series B}
  \bold{66} (2004) 79-83.
}
\seealso{
  Mark correlation \code{\link{markcorr}},
  mark variogram \code{\link{markvario}} for numeric marks.
  
  Mark connection function \code{\link{markconnect}} and 
  multitype K-functions \code{\link{Kcross}}, \code{\link{Kdot}}
  for factor-valued marks.
}
\examples{
    plot(Emark(spruces))
    E <- Emark(spruces, method="density", kernel="epanechnikov")
    plot(Vmark(spruces))
}
\author{
  Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  and
  Rolf Turner \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{nonparametric}


