\name{Smooth}
\alias{Smooth}
\title{Spatial smoothing of data}
\description{
  Generic function to perform spatial smoothing of spatial data.
}
\usage{
  Smooth(X, ...)
}
\arguments{
  \item{X}{Some kind of spatial data}
  \item{\dots}{Arguments passed to methods.}
}
\details{
  This generic function calls an appropriate method
  to perform spatial smoothing on the spatial dataset \code{X}.

  Methods for this function include
  \itemize{
    \item \code{\link{Smooth.ppp}} for point patterns
    \item \code{\link{Smooth.msr}} for measures
    \item \code{\link{Smooth.fv}} for function value tables
  }
}
\seealso{
  \code{\link{Smooth.ppp}},
  \code{\link{Smooth.im}},
  \code{\link{Smooth.msr}},
  \code{\link{Smooth.fv}}.
}
\value{
  An object containing smoothed values of the input data,
  in an appropriate format. See the documentation for the methods.
}
\author{
  \adrian
  and \rolf
}
\keyword{spatial}
\keyword{methods}
\keyword{smooth}
