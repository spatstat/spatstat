\name{markmarkscatter}
\alias{markmarkscatter}
\title{
  Mark-Mark Scatter Plot
}
\description{
  Generates the mark-mark scatter plot of a point pattern.
}
\usage{
markmarkscatter(X, rmax, \dots, col = NULL, symap = NULL, transform=I, jit=FALSE)
}
\arguments{
  \item{X}{
    A point pattern (object of class \code{"ppp"}, \code{"pp3"},
    \code{"lpp"} or \code{"ppx"}) with numeric marks.
  }
  \item{rmax}{
    Maximum distance between pairs of points
    which contribute to the plot.
  }
  \item{\dots}{
    Additional arguments passed to \code{\link{plot.ppp}} to
    control the scatterplot.
  }
  \item{transform}{
    Optional. A function which should be applied to the mark values.
  }
  \item{jit}{
    Logical value indicating whether mark values should be randomly
    perturbed using \code{\link[base]{jitter}}. 
  }
  \item{col}{
    Optional. A vector of colour values, or a \code{\link{colourmap}}
    to be used to portray the pairwise distance values. Ignored
    if \code{symap} is given.
  }
  \item{symap}{
    Optional. A \code{\link{symbolmap}} to be used to
    portray the pairwise distance values. Overrides \code{col}.
  }
}
\details{
  The mark-mark scatter plot (Ballani et al, 2019) is a
  scatterplot of the mark values of all pairs of distinct points in \code{X}
  which are closer than the distance \code{rmax}.
  The dots in the scatterplot are coloured according to the
  pairwise distance between the two spatial points.
  The plot is augmented by three curves explained by
  Ballani et al (2019).

  If the marks only take a few different values, then
  it is usually appropriate to apply random perturbation (jitter)
  to the mark values, by setting \code{jit=TRUE}.
}
\value{
  Null.
}
\references{
  Ballani, F., Pommerening, A. and Stoyan, D. (2019)
  Mark-mark scatterplots improve pattern analysis in spatial plant
  ecology.
  \emph{Ecological Informatics} \bold{49}, 13--21.
}
\author{
  Adrian Baddeley (coded from the description in Ballani et al.)
}
\examples{
  markmarkscatter(longleaf, 10)

  markmarkscatter(spruces, 10, jit=TRUE)
}
\keyword{spatial}
\keyword{hplot}
