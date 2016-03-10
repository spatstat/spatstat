\name{nnclean}
\alias{nnclean}
\alias{nnclean.ppp}
\alias{nnclean.pp3}
\title{
  Nearest Neighbour Clutter Removal
}
\description{
  Detect features in a 2D or 3D spatial point pattern
  using nearest neighbour clutter removal.
}
\usage{
  nnclean(X, k, ...)

  \method{nnclean}{ppp}(X, k, ..., 
                     edge.correct = FALSE, wrap = 0.1,
                     convergence = 0.001, plothist = FALSE,
                     verbose = TRUE, maxit = 50)

  \method{nnclean}{pp3}(X, k, ..., 
                     convergence = 0.001, plothist = FALSE,
                     verbose = TRUE, maxit = 50)
}
\arguments{
  \item{X}{
    A two-dimensional spatial point pattern (object of class
    \code{"ppp"}) or a three-dimensional point pattern
    (object of class \code{"pp3"}).
  }
  \item{k}{
    Degree of neighbour: \code{k=1} means nearest neighbour,
    \code{k=2} means second nearest, etc.
  }
  \item{\dots}{
    Arguments passed to \code{\link{hist.default}} to control
    the appearance of the histogram, if \code{plothist=TRUE}.
  }
  \item{edge.correct}{
    Logical flag specifying whether periodic edge correction
    should be performed (only implemented in 2 dimensions).
  }
  \item{wrap}{
    Numeric value specifying the relative size of the margin
    in which data will be replicated for the 
    periodic edge correction (if \code{edge.correct=TRUE}).
    A fraction of window width and window height.
  }
  \item{convergence}{
    Relative tolerance threshold for testing convergence of EM algorithm.
  }
  \item{maxit}{
    Maximum number of iterations for EM algorithm.
  }
  \item{plothist}{
    Logical flag specifying whether to plot a diagnostic histogram
    of the nearest neighbour distances and the fitted distribution.
  }
  \item{verbose}{
    Logical flag specifying whether to print progress reports.
  }
}
\details{
  Byers and Raftery (1998) developed a technique for recognising
  features in a spatial point pattern in the presence of
  random clutter.

  For each point in the pattern, the distance to the
  \eqn{k}th nearest neighbour is computed. Then the E-M algorithm is
  used to fit a mixture distribution to the
  \eqn{k}th nearest neighbour distances.
  The mixture components represent the feature and the clutter. The
  mixture model can be used to classify each point as belong to one
  or other component.

  The function \code{nnclean} is generic, with methods for
  two-dimensional point patterns (class \code{"ppp"})
  and three-dimensional point patterns (class \code{"pp3"})
  currently implemented.

  The result is a point pattern (2D or 3D) with two additional
  columns of marks:
  \describe{
    \item{class}{
      A factor, with levels \code{"noise"} and \code{"feature"},
      indicating the maximum likelihood classification of each point.
    }
    \item{prob}{
      Numeric vector giving the estimated probabilities
      that each point belongs to a feature.
    }
  }
  The object also has extra information stored in attributes:
  \code{"theta"} contains the fitted parameters
  of the mixture model, \code{"info"} contains
  information about the fitting procedure, and \code{"hist"} contains
  the histogram structure returned from \code{\link{hist.default}}
  if \code{plothist = TRUE}.
}
\value{
  An object of the same kind as \code{X},
  obtained by attaching marks to the points of \code{X}.

  The object also has attributes, as described under Details.
}
\references{
  Byers, S. and Raftery, A.E. (1998)
  Nearest-neighbour clutter removal for estimating features
  in spatial point processes.
  \emph{Journal of the American Statistical Association}
  \bold{93}, 577--584.
}
\author{
  Original by Simon Byers and Adrian Raftery.
  Adapted for \pkg{spatstat} by Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}.
}
\seealso{
  \code{\link{nndist}}, 
  \code{\link{split.ppp}}, 
  \code{\link{cut.ppp}}
}
\examples{
  data(shapley)
  X <- nnclean(shapley, k=17)
  plot(X, which.marks=1, chars=c(".", "+"), cols=1:2)
  plot(X, which.marks=2, cols=function(x)hsv(0.2+0.8*(1-x),1,1))
  Y <- split(X, un=TRUE)
  plot(Y, chars="+", cex=0.5)
  marks(X) <- marks(X)$prob
  plot(cut(X, breaks=3), chars=c(".", "+", "+"), cols=1:3)
}
\keyword{spatial}
\keyword{classif}
