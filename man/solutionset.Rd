\name{solutionset}
\alias{solutionset}
\title{Evaluate Logical Expression Involving Pixel Images and Return
  Region Where Expression is True}
\description{
  Given a logical expression involving one or more pixel images,
  find all pixels where the expression is true,
  and assemble these pixels into a window.
}
\usage{
  solutionset(\dots, envir)
}
\arguments{
  \item{\dots}{An expression in the \R language, involving one or more
    pixel images.}
  \item{envir}{Optional. The environment in which to evaluate the
    expression.}
}
\details{
  Given a logical expression involving one or more pixel images,
  this function will find all pixels where the expression is true,
  and assemble these pixels into a spatial window.

  Pixel images in \code{spatstat}
  are represented by objects of class \code{"im"}
  (see \code{\link{im.object}}). These are essentially matrices of
  pixel values, with extra attributes recording the pixel dimensions,
  etc.

  Suppose \code{X} is a pixel image. Then \code{eval.im(abs(X) > 3)}
  will find all the pixels in \code{X} for which the pixel value
  is greater than 3 in absolute value, and return a window containing
  all these pixels.

  Suppose \code{X} and \code{Y} are two pixel images with compatible
  dimensions: they have the same number of pixels, the same physical
  size of pixels, and the same bounding box. Then
  \code{eval.im(X > Y)} will find all pixels for which the
  pixel value of \code{X} is greater than the corresponding pixel value
  of \code{Y}, and return a window containing these pixels.
  
  In general, \code{expr} can be any logical expression involving
  (a) the \emph{names} of pixel images, (b) scalar
  constants, and (c) functions which are vectorised.
  See the Examples.

  The expression \code{expr} is evaluated by \code{\link{eval.im}}.
  The expression \code{expr} must be vectorised.
  There must be at least one pixel image in the expression.
  All images must have compatible dimensions.
}
\value{
  A spatial window
  (object of class \code{"owin"}, see \code{\link{owin.object}}).
}
\seealso{
  \code{\link{im.object}},
  \code{\link{owin.object}},
  \code{\link{eval.im}},
  \code{\link{levelset}}
}
\examples{
  # test images
  X <- as.im(function(x,y) { x^2 - y^2 }, unit.square())
  Y <- as.im(function(x,y) { 3 * x + y  - 1}, unit.square())

  W <- solutionset(abs(X) > 0.1)
  W <- solutionset(X > Y)
  W <- solutionset(X + Y >= 1)

  area.owin(solutionset(X < Y))
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{programming}
\keyword{manip}
