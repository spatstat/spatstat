\name{eval.fasp}
\alias{eval.fasp}
\title{Evaluate Expression Involving Function Arrays}
\description{
  Evaluates any expression involving one or more function arrays
  (\code{fasp} objects) and returns another function array.
}
\usage{
  eval.fasp(expr, envir, dotonly=TRUE)
}
\arguments{
  \item{expr}{
    An expression involving the names of objects of class \code{"fasp"}.
  }
  \item{envir}{
    Optional. The environment in which to evaluate the expression.
  }
  \item{dotonly}{Logical. Passed to \code{\link{eval.fv}}.}
}
\details{
  This is a wrapper to make it easier to perform
  pointwise calculations with the arrays of summary functions
  used in spatial statistics.

  A function array (object of class \code{"fasp"}) can be regarded as a matrix
  whose entries are functions. Objects of this kind
  are returned by the command \code{\link{alltypes}}.

  Suppose \code{X} is an object of class \code{"fasp"}.
  Then \code{eval.fasp(X+3)} effectively adds 3 to the value of
  every function in the array \code{X}, and returns
  the resulting object. 

  Suppose \code{X} and \code{Y} are two objects of class \code{"fasp"}
  which are compatible (for example the arrays
  must have the same dimensions). Then 
  \code{eval.fasp(X + Y)} will add the corresponding functions in
  each cell of the arrays \code{X} and \code{Y},
  and return the resulting array of functions.

  Suppose \code{X} is an object of class \code{"fasp"}
  and \code{f} is an object of class \code{"fv"}.
  Then \code{eval.fasp(X + f)} will add the function \code{f}
  to the functions in each cell of the array \code{X},
  and return the resulting array of functions.

  In general, \code{expr} can be any expression involving
  (a) the \emph{names} of objects of class \code{"fasp"} or \code{"fv"},
  (b) scalar constants, and (c) functions which are vectorised.
  See the Examples.

  First \code{eval.fasp} determines which of the \emph{variable names}
  in the expression \code{expr} refer to objects of class \code{"fasp"}.
  The expression is then evaluated for each cell of the array
  using \code{\link{eval.fv}}.

  The expression \code{expr} must be vectorised.
  There must be at least one object of class \code{"fasp"} in the expression.
  All such objects must be compatible.
}
\value{
  Another object of class \code{"fasp"}.
}
\seealso{
  \code{\link{fasp.object}},
  \code{\link{Kest}}
}
\examples{
  # manipulating the K function
  K <- alltypes(amacrine, "K")

  # expressions involving a fasp object
  eval.fasp(K + 3)
  L <- eval.fasp(sqrt(K/pi))

  # expression involving two fasp objects
  D <- eval.fasp(K - L)

  # subtracting the unmarked K function from the cross-type K functions
  K0 <- Kest(unmark(amacrine))
  DK <- eval.fasp(K - K0)
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{manip}
\keyword{programming}
