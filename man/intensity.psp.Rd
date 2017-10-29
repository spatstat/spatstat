\name{intensity.psp}  
\alias{intensity.psp}
\title{
  Empirical Intensity of Line Segment Pattern
}
\description{
  Computes the average total length of segments per unit area
  in a spatial pattern of line segments.
}
\usage{
\method{intensity}{psp}(X, ..., weights=NULL)
}
\arguments{
  \item{X}{
    A line segment pattern (object of class \code{"psp"}).
  }
  \item{weights}{
    Optional.
    Numeric vector of weights attached to the segments of \code{X}.
    Alternatively, an \code{expression} which can be evaluated to
    give a vector of weights.
  }
  \item{\dots}{
    Ignored.
  }
}
\details{
  This is a method for the generic function \code{\link{intensity}}.
  It computes the empirical intensity of a line segment pattern
  (object of class \code{"psp"}),
  i.e. the average total segment length per unit area.

  If the segment pattern is multitype, the intensities of the
  different types are computed separately.

  Note that the intensity will be computed as the length per area
  in units per square unit, based on the unit of length for \code{X},
  given by \code{unitname(X)}. If the unit of length is a strange multiple
  of a standard unit, like \code{5.7 metres}, then it can be converted
  to the standard unit using \code{\link{rescale}}. See the Examples.

  If \code{weights} are given, then the intensity is
  computed as the total \emph{weight times length} per square unit.
  The argument \code{weights} should be a numeric vector
  of weights for each point of \code{X} (weights may be negative or
  zero). 

  Alternatively \code{weights} can be an \code{expression}
  which will be evaluated for the dataset to yield a vector of weights.
  The expression may involve the Cartesian coordinates \eqn{x,y} of the
  points, and the marks of the points, if any. Variable names
  permitted in the expression include \code{x0}, \code{x1}, \code{y0},
  \code{y1} for the coordinates of the segment endpoint,
  the name \code{marks} if \code{X} has a single column of marks,
  the names of any columns of marks if \code{X} has a data frame of
  marks, and the names of constants or functions that exist
  in the global environment. See the Examples.
}
\value{
  A numeric value (giving the intensity) or numeric vector
  (giving the intensity for each possible type).
}
\seealso{
  \code{\link{intensity}}
}
\examples{
  S <- as.psp(simplenet)
  intensity(S)
  intensity(S, weights=runif(nsegments(S)))
  intensity(S, weights=expression((x0+x1)/2))
}
\author{
  \spatstatAuthors.
}
\keyword{spatial}
\keyword{nonparametric}
