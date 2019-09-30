\name{bw.voronoi}
\alias{bw.voronoi}
\title{
  Cross Validated Bandwidth Selection for Voronoi Estimator of Intensity
  on a Network
}
\description{
  Uses cross-validation to select a smoothing bandwidth
  for the Voronoi estimate of point process intensity
  on a linear network.
}
\usage{
bw.voronoi(X, \dots, probrange = c(0.2, 0.8), nprob = 10,
           prob = NULL, nrep = 100, verbose = TRUE, warn=TRUE)
}
\arguments{
  \item{X}{
    Point pattern on a linear network (object of class \code{"lpp"}).
  }
  \item{\dots}{
    Ignored.
  }
  \item{probrange}{
    Numeric vector of length 2 giving the range of bandwidths
    (retention probabilities)
    to be assessed.
  }
  \item{nprob}{
    Integer. Number of bandwidths to be assessed.
  }
  \item{prob}{
    Optional. A numeric vector of bandwidths (retention probabilities)
    to be assessed. Entries must be probabilities between 0 and 1.
    Overrides \code{nprob} and \code{probrange}.
  }
  \item{nrep}{
    Number of simulated realisations to be used for the computation.
  }
  \item{verbose}{
    Logical value indicating whether to print progress reports.
  }
  \item{warn}{
    Logical. If \code{TRUE}, issue a warning if the maximum of
    the cross-validation criterion occurs at one of the ends of the
    search interval.
  }
}
\details{
  This function uses likelihood cross-validation to choose the optimal value of the
  thinning fraction \code{f} (the retention probability)
  to be used in the smoothed Voronoi estimator of point process
  intensity \code{\link{densityVoronoi.lpp}}.
}
\value{
  A numerical value giving the selected bandwidth.
  The result also belongs to the class \code{"bw.optim"}
  which can be plotted.
}
\references{
  Moradi, M., Cronie, 0., Rubak, E., Lachieze-Rey, R.,
  Mateu, J. and Baddeley, A. (2019)
  Resample-smoothing of Voronoi intensity estimators.
  \emph{Statistics and Computing}, in press.
}
\author{
  \spatstatAuthors and Mehdi Moradi.
}
\seealso{
  \code{\link{densityVoronoi.lpp}}
}
\examples{
   np <- if(interactive()) 10 else 3
   nr <- if(interactive()) 100 else 2
   b <- bw.voronoi(spiders, nprob=np, nrep=nr)
   b
   plot(b)
}
\keyword{spatial}
\keyword{methods}
\keyword{smooth}
