\name{eval.linim}
\alias{eval.linim}
\title{Evaluate Expression Involving Pixel Images on Linear Network}
\description{
  Evaluates any expression involving one or more pixel images
  on a linear network, and returns a pixel image on the same linear network.
}
\usage{
  eval.linim(expr, envir, harmonize=TRUE)
}
\arguments{
  \item{expr}{An expression in the \R language,
    involving the names of objects of class \code{"linim"}.}
  \item{envir}{Optional. The environment in which to evaluate the
    expression.}
  \item{harmonize}{
    Logical. Whether to resolve inconsistencies between
    the pixel grids.
  }
}
\details{
  This function a wrapper to make it easier to perform
  pixel-by-pixel calculations. It is one of several functions
  whose names begin with \code{eval} which work on objects of
  different types. This particular function is designed to work with
  objects of class \code{"linim"} which represent
  pixel images on a linear network. 

  Suppose \code{X} is a pixel image on a linear network (object of
  class \code{"linim"}. Then \code{eval.linim(X+3)}
  will add 3 to the value of every pixel in \code{X}, and return
  the resulting pixel image on the same linear network.

  Suppose \code{X} and \code{Y} are two pixel images on the same
  linear network, with compatible pixel dimensions. Then
  \code{eval.linim(X + Y)} will add the corresponding pixel values in
  \code{X} and \code{Y}, and return the resulting pixel image
  on the same linear network.

  In general, \code{expr} can be any expression in the R language involving
  (a) the \emph{names} of pixel images, (b) scalar
  constants, and (c) functions which are vectorised.
  See the Examples.

  First \code{eval.linim} determines which of the \emph{variable names}
  in the expression \code{expr} refer to pixel images. Each such name
  is replaced by a matrix containing the pixel values. The expression is
  then evaluated. The result should be a matrix; it is taken as
  the matrix of pixel values.
  
  The expression \code{expr} must be vectorised.
  There must be at least one linear pixel image in the expression.

  All images must have compatible dimensions.
  If \code{harmonize=TRUE}, images that have incompatible dimensions
  will be resampled so that they are compatible.
  If \code{harmonize=FALSE}, images that are incompatible will cause an error.
}
\value{
  An image object of class \code{"linim"}.
}
\seealso{
  \code{\link{eval.im}},
  \code{\link{linim}}
}
\examples{
  example(linim)

  Y <- linim(letterA, as.im(function(x,y){y^2+x}, W=M))
  
  eval.linim(X + 3)
  eval.linim(X - Y)
  eval.linim(abs(X - Y))
  Z <- eval.linim(sin(X * pi) + Y)
}
\author{Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  and Rolf Turner \email{r.turner@auckland.ac.nz}
  
}
\keyword{spatial}
\keyword{manip}
\keyword{programming}
