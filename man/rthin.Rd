\name{rthin}
\alias{rthin}
\title{Random Thinning}
\description{
  Applies independent random thinning to a point pattern.
}
\usage{
rthin(X, P, \dots, nsim=1, drop=TRUE)
}
\arguments{
  \item{X}{
    A point pattern (object of class \code{"ppp"})
    that will be thinned.
  }
  \item{P}{
    Data giving the retention probabilities, i.e. the probability
    that each point in \code{X} will be retained.
    Either a single number, or a vector of numbers,
    or a \code{function(x,y)}, or a pixel image (object of class
    \code{"im"}).
  }
  \item{\dots}{
    Additional arguments passed to \code{P}, if it is a function.
  }
  \item{nsim}{Number of simulated realisations to be generated.}
  \item{drop}{
    Logical. If \code{nsim=1} and \code{drop=TRUE} (the default), the
    result will be a point pattern, rather than a list 
    containing a point pattern.
  }
}
\value{
  A point pattern (an object of class \code{"ppp"})
  if \code{nsim=1}, or a list of point patterns if \code{nsim > 1}.
}
\details{
  In a random thinning operation, each point of the pattern \code{X}
  is randomly either deleted or retained (i.e. not deleted).
  The result is a point pattern,
  consisting of those points of \code{X} that were retained.

  Independent random thinning means that the retention/deletion of each
  point is independent of other points.

  The argument \code{P} determines the probability of \bold{retaining}
  each point. It may be
  \describe{
    \item{a single number,}{so that each point will be retained with
      the same probability \code{P};
    }
    \item{a vector of numbers,}{so that the \code{i}th point of \code{X}
      will be retained with probability \code{P[i]};
    }
    \item{a function \code{P(x,y)},}{so that a point at a location
      \code{(x,y)} will be retained with probability \code{P(x,y)};
    }
    \item{a pixel image,}{containing values of the retention probability
      for all locations in a region encompassing the point pattern.
    }
  }
  If \code{P} is a function, it should be \sQuote{vectorised}, that is,
  it should accept vector arguments \code{x,y} and should yield a
  numeric vector of the same length. The function may have extra
  arguments which are passed through the \code{\dots} argument.
}
\section{Reproducibility}{
  The algorithm for random thinning was changed in \pkg{spatstat}
  version \code{1.42-3}. Set \code{spatstat.options(fastthin=FALSE)}
  to use the previous, slower algorithm, if it is desired to reproduce
  results obtained with earlier versions.
}
\examples{
  data(redwood)
  plot(redwood, main="thinning")
  
  # delete 20\% of points
  Y <- rthin(redwood, 0.8)
  points(Y, col="green", cex=1.4)

  # function
  f <- function(x,y) { ifelse(x < 0.4, 1, 0.5) }
  Y <- rthin(redwood, f)

  # pixel image
  Z <- as.im(f, Window(redwood))
  Y <- rthin(redwood, Z)
}
\author{\adrian
  
  
  and \rolf
  
}
\keyword{spatial}
\keyword{datagen}
