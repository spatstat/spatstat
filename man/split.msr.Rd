\name{split.msr}
\alias{split.msr}
\title{
  Divide a Measure into Parts
}
\description{
  Decomposes a measure into components, each component being a measure.
}
\usage{
  \method{split}{msr}(x, f, drop = FALSE, \dots)
}
\arguments{
  \item{x}{
    Measure (object of class \code{"msr"}) to be decomposed.
  }
  \item{f}{
    Factor or tessellation determining the decomposition.
    Argument passed to \code{\link{split.ppp}}.
    See Details.
  }
  \item{drop}{
    Logical value indicating whether empty components should be retained
    in the list (\code{drop=FALSE}, the default) or deleted (\code{drop=TRUE}).
  }
  \item{\dots}{
    Ignored.
  }
}
\details{
  An object of class \code{"msr"} represents a signed (i.e. real-valued)
  or vector-valued measure in the \pkg{spatstat} package.
  See \code{\link{msr}} for explanation.

  This function is a method for the generic
  \code{\link[base]{split}}. It divides the measure \code{x} into
  components, each of which is a measure.

  A measure \code{x} is represented in \pkg{spatstat} by 
  a finite set of sample points with values attached to them.
  The function \code{split.msr} divides this pattern of sample points
  into several sub-patterns of points using \code{\link{split.ppp}}.
  For each sub-pattern, the values attached to these points are
  extracted from \code{x}, and these values and sample points
  determine a measure, which is a
  component or piece of the original \code{x}.

  The argument \code{f} can be missing, if the sample points of \code{x}
  are multitype points. In this case, \code{x} represents a measure
  associated with marked spatial locations, and the command \code{split(x)}
  separates \code{x} into a list of component measures, one for each
  possible mark.

  Otherwise the argument \code{f} is passed to \code{\link{split.ppp}}.
  It should be either a factor (of length equal to the number of
  sample points of \code{x}) or a tessellation (object of class
  \code{"tess"} representing a division of space into tiles)
  as documented under \code{\link{split.ppp}}.
}
\value{
  A list, each of whose entries is a measure (object of class
  \code{"msr"}). 
}
\author{
  \spatstatAuthors
}
\seealso{
  \code{\link{msr}},
  \code{\link{[.msr}},
  \code{\link{with.msr}}
}
\examples{
  ## split by tessellation
  a <- residuals(ppm(cells ~ x))
  aa <- split(a, dirichlet(runifpoint(4)))
  aa
  sapply(aa, integral)

  ## split by type of point
  b <- residuals(ppm(amacrine ~ marks + x))
  bb <- split(b)
  bb
}
\keyword{spatial}
\keyword{manip}
