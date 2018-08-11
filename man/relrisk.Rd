\name{relrisk}
\alias{relrisk}
\title{
  Estimate of Spatially-Varying Relative Risk
}
\description{
  Generic command to estimate the 
  spatially-varying probability of each type of point, or the ratios of
  such probabilities.
}
\usage{
  relrisk(X, \dots) 
}
\arguments{
  \item{X}{
    Either a point pattern (class \code{"ppp"})
    or a fitted point process model (class \code{"ppm"})
    from which the probabilities will be estimated.
  }
  \item{\dots}{
    Additional arguments appropriate to the method.
  }
}
\details{
  In a point pattern containing several different types of points,
  we may be interested in the spatially-varying probability of each
  possible type, or the relative risks which are the
  ratios of such probabilities.
  
  The command \code{\link{relrisk}} is generic and can be used to
  estimate relative risk in different ways.
  
  The function \code{\link{relrisk.ppp}} is the method for point pattern
  datasets. It computes \emph{nonparametric} estimates of relative risk
  by kernel smoothing.

  The function \code{\link{relrisk.ppm}} is the method for fitted point
  process models (class \code{"ppm"}). It computes \emph{parametric}
  estimates of relative risk, using the fitted model.
}
\seealso{
  \code{\link{relrisk.ppp}},
  \code{\link{relrisk.ppm}}.
}
\author{\adrian
  
  
  \rolf
  
  and \ege
  
}
\keyword{spatial}
