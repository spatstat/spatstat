\name{quadratresample}
\alias{quadratresample}
\title{Resample a Point Pattern by Resampling Quadrats}
\description{
  Given a point pattern dataset, create a resampled point pattern
  by dividing the window into rectangular quadrats and randomly
  resampling the list of quadrats. 
}
\usage{
quadratresample(X, nx, ny=nx, ...,
                replace = FALSE, nsamples = 1,
                verbose = (nsamples > 1))
}
\arguments{
  \item{X}{
    A point pattern dataset (object of class \code{"ppp"}).
  }
  \item{nx,ny}{
    Numbers of quadrats in the \eqn{x} and \eqn{y} directions.
  }
  \item{\dots}{Ignored.}
  \item{replace}{
    Logical value. Specifies whether quadrats should be sampled
    with or without replacement.
  }
  \item{nsamples}{Number of randomised point patterns to be generated.}
  \item{verbose}{Logical value indicating whether to print progress reports.}
}
\details{
  This command implements a very simple bootstrap resampling procedure
  for spatial point patterns \code{X}.

  The dataset \code{X} must be a point pattern (object of class
  \code{"ppp"}) and its observation window must be a rectangle.

  The window is first divided into \code{N = nx * ny} rectangular tiles
  (quadrats) of equal size and shape. 
  To generate one resampled point pattern, a random sample of
  \code{N} quadrats is selected from the list of \code{N} quadrats,
  with replacement (if \code{replace=TRUE}) or without replacement
  (if \code{replace=FALSE}). The \eqn{i}th quadrat in the original
  dataset is then replaced by the \eqn{i}th sampled quadrat, after the
  latter is shifted so that it
  occupies the correct spatial position. The quadrats are then
  reconstituted into a point pattern inside the same window as \code{X}.

  If \code{replace=FALSE}, this procedure effectively involves a random
  permutation of the quadrats. The resulting resampled point pattern has
  the same number of points as \code{X}.
  If \code{replace=TRUE}, the number of points in the resampled point
  pattern is random.
}
\value{
  A point pattern (if \code{nsamples = 1}) or a 
  list of point patterns (if \code{nsamples > 1}).
}
\author{\adrian
  
  
  and \rolf
  
}
\seealso{
  \code{\link{quadrats}},
  \code{\link{quadratcount}}.

  See \code{\link{varblock}} to estimate the variance of
  a summary statistic by block resampling.
}
\examples{
  data(bei)
  quadratresample(bei, 6, 3)
}
\keyword{spatial}
\keyword{datagen}
