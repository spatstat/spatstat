\name{with.ssf}
\alias{with.ssf}
\alias{apply.ssf}
\title{
  Evaluate Expression in a Spatially Sampled Function
}
\description{
  Given a spatially sampled function,
  evaluate an expression involving the function values.
}
\usage{
  apply.ssf(X, \dots)

  \method{with}{ssf}(data, \dots)
}
\arguments{
  \item{X, data}{
    A spatially sampled function (object of class \code{"ssf"}).
  }
  \item{\dots}{
    Arguments passed to \code{\link{with.default}} or
    \code{\link{apply}} specifying what to compute.
  }
}
\details{
  An object of class \code{"ssf"} represents a
  function (real- or vector-valued) that has been
  sampled at a finite set of points.
  It contains a data frame
  which provides the function values
  at the sample points. 

  In \code{with.ssf}, the expression specified by \code{\dots}
  will be evaluated in this dataframe.
  In \code{apply.ssf}, the dataframe will be subjected to
  the \code{\link{apply}} operator using the additional arguments
  \code{\dots}. 
  
  If the result of evaluation
  is a data frame with one row for each data point,
  or a numeric vector with one entry for each data point,
  then the result will be an object of class \code{"ssf"}
  containing this information. Otherwise, the result will be
  a numeric vector.
}
\value{
  An object of class \code{"ssf"} or a numeric vector.
}
\author{
  \adrian.
}
\seealso{
 \code{\link{ssf}}
}
\examples{
  a <- ssf(cells, data.frame(d=nndist(cells), i=1:npoints(cells)))
  with(a, i/d)
}
\keyword{spatial}
\keyword{manip}
\keyword{programming}

