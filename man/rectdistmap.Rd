\name{rectdistmap}
\alias{rectdistmap}
\title{
  Distance Map Using Rectangular Distance Metric
}
\description{
  Computes the distance map of a spatial region
  based on the rectangular distance metric.
}
\usage{
rectdistmap(X, asp = 1, npasses=1, verbose=FALSE)
}
\arguments{
  \item{X}{
    A window (object of class \code{"owin"}).
  }
  \item{asp}{
    Aspect ratio for the metric.
    See Details.
  }
  \item{npasses}{
    Experimental. 
  }
  \item{verbose}{
    Logical value indicating whether to print trace information.
  }
}
\details{
  This function computes the distance map of the
  spatial region \code{X} using the rectangular distance metric
  with aspect ratio \code{asp}. This metric is defined so that
  the set of all points lying at most 1 unit away from the origin
  (according to the metric) form a rectangle of width 1
  and height \code{asp}.
}
\value{
  A pixel image (object of class \code{"im"}).
}
\author{
  \adrian
}
\seealso{
  \code{\link{distmap}}
}
\examples{
  V <- letterR
  Frame(V) <- grow.rectangle(Frame(V), 0.5)
  plot(rectdistmap(V))
}
\keyword{spatial}
\keyword{math}
