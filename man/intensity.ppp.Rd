\name{intensity.ppp}  %DontDeclareMethods
\alias{intensity.ppp}
\alias{intensity.splitppp}
\title{
  Empirical Intensity of Point Pattern
}
\description{
  Computes the average number of points per unit area
  in a point pattern dataset.
}
\usage{
\method{intensity}{ppp}(X, ..., weights=NULL)

\method{intensity}{splitppp}(X, ..., weights=NULL)
}
\arguments{
  \item{X}{
    A point pattern (object of class \code{"ppp"}).
  }
  \item{weights}{
    Optional.
    Numeric vector of weights attached to the points of \code{X}.
    Alternatively, an \code{expression} which can be evaluated to
    give a vector of weights.
  }
  \item{\dots}{
    Ignored.
  }
}
\details{
  This is a method for the generic function \code{\link{intensity}}.
  It computes the empirical intensity of a point pattern
  (object of class \code{"ppp"}),
  i.e. the average density of points per unit area.

  If the point pattern is multitype, the intensities of the
  different types are computed separately.

  Note that the intensity will be computed as the number of points
  per square unit, based on the unit of length for \code{X},
  given by \code{unitname(X)}. If the unit of length is a strange multiple
  of a standard unit, like \code{5.7 metres}, then it can be converted
  to the standard unit using \code{\link{rescale}}. See the Examples.

  If \code{weights} are given, then the intensity is
  computed as the total \emph{weight} per square unit.
  The argument \code{weights} should be a numeric vector
  of weights for each point of \code{X} (weights may be negative or
  zero).

  Alternatively \code{weights} can be an \code{expression}
  which will be evaluated for the dataset to yield a vector of weights.
  The expression may involve the Cartesian coordinates \eqn{x,y} of the
  points, and the marks of the points, if any. Variable names
  permitted in the expression include \code{x} and \code{y},
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
  \code{\link{intensity}}, 
  \code{\link{intensity.ppm}}
}
\examples{
  japanesepines
  intensity(japanesepines)
  unitname(japanesepines)
  intensity(rescale(japanesepines))

  intensity(amacrine)
  intensity(split(amacrine))

  # numeric vector of weights
  volumes <- with(marks(finpines), (pi/4) * height * diameter^2)
  intensity(finpines, weights=volumes)

  # expression for weights
  intensity(finpines, weights=expression((pi/4) * height * diameter^2))
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{nonparametric}
