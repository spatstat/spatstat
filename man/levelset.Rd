\name{levelset}
\alias{levelset}
\title{Level Set of a Pixel Image}
\description{
  Given a pixel image, find all pixels which have values less than a
  specified threshold value (or greater than a threshold, etc),
  and assemble these pixels into a window.
}
\usage{
  levelset(X, thresh, compare="<=")
}
\arguments{
  \item{X}{A pixel image (object of class "im")}.
  \item{thresh}{Threshold value.
    A single number or value compatible with the pixel values in
    \code{X}}.
  \item{compare}{Character string specifying one of the comparison
    operators \code{"<", ">", "==", "<=", ">=", "!="}. 
  }
}
\details{
  If \code{X} is a pixel image with numeric values,
  then \code{levelset(X, thresh)} finds the region of space
  where the pixel values are less than or equal to
  the threshold value \code{thresh}.
  This region is returned as a spatial window.

  The argument \code{compare} specifies how the pixel values should be
  compared with the threshold value.
  Instead of requiring pixel values to be less than or equal to
  \code{thresh}, you can specify that they must be less than (\code{<}),
  greater than (\code{>}), equal to (\code{==}), greater than or equal
  to (\code{>=}), or not equal to (\code{!=}) the threshold value
  \code{thresh}.

  If \code{X} has non-numeric pixel values (for example, logical
  or factor values) it is advisable to use only the comparisons
  \code{==} and \code{!=}, unless you really know what you are doing.

  For more complicated logical comparisons, see \code{\link{solutionset}}.
}
\value{
  A spatial window (object of class \code{"owin"},
  see \code{\link{owin.object}}) containing the pixels satisfying the
  constraint.
}
\seealso{
  \code{\link{im.object}},
  \code{\link{as.owin}},
  \code{\link{solutionset}}.
}
\examples{
  # test image
  X <- as.im(function(x,y) { x^2 - y^2 }, unit.square())

  W <- levelset(X, 0.2)
  W <- levelset(X, -0.3, ">")

  # compute area of level set
  area.owin(levelset(X, 0.1))
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{programming}
\keyword{manip}
