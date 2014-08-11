\name{im}
\alias{im}
\title{Create a Pixel Image Object}
\description{
  Creates an object of
  class \code{"im"} representing a two-dimensional pixel image.
}
\usage{
  im(mat, xcol=seq_len(ncol(mat)), yrow=seq_len(nrow(mat)),
   xrange=NULL, yrange=NULL,
   unitname=NULL)
}
\arguments{
  \item{mat}{
    matrix or vector containing the pixel values of the image.
  }
  \item{xcol}{
    vector of \eqn{x} coordinates for the pixel gid
  }
  \item{yrow}{
    vector of \eqn{y} coordinates for the pixel grid
  }
  \item{xrange,yrange}{
    Optional. Vectors of length 2 giving the \eqn{x} and \eqn{y}
    limits of the enclosing rectangle.
    (Ignored if \code{xcol}, \code{yrow} are present.)
  }
  \item{unitname}{
    Optional. Name of unit of length. Either a single character string,
    or a vector of two character strings giving the
    singular and plural forms, respectively.
  }
}
\details{
  This function creates an object of class \code{"im"} representing
  a two-dimensional pixel image. See \code{\link{im.object}}
  for details of this class.

  The matrix \code{mat} contains the \sQuote{greyscale} values
  for a rectangular grid of pixels.
  Note carefully that the entry \code{mat[i,j]}
  gives the pixel value at the location \code{(xcol[j],yrow[i])}.
  That is, the \bold{row} index of the matrix \code{mat} corresponds
  to increasing \bold{y} coordinate, while the column index of \code{mat}
  corresponds to increasing \bold{x} coordinate.
  Thus \code{yrow} has one entry for each row of \code{mat}
  and \code{xcol} has one entry for each column of \code{mat}.
  Under the usual convention in \R, a correct
  display of the image would be obtained by transposing the matrix, e.g.
  \code{image.default(xcol, yrow, t(mat))}, if you wanted to do it by hand.

  The entries of \code{mat} may be numeric (real or integer), complex, 
  logical, character, or factor values.
  If \code{mat} is not a matrix, it will be converted into
  a matrix with \code{nrow(mat) = length(yrow)} and
  \code{ncol(mat) = length(xcol)}.
  
  To make a factor-valued image, note that 
  \R has a quirky way of handling matrices with
  factor-valued entries. The command \code{\link{matrix}} cannot be used
  directly, because it destroys factor information.
  To make a factor-valued image, do one of the following:
  \itemize{
    \item
    Create a \code{factor} containing the pixel values,
    say \code{mat <- factor(.....)}, 
    and then assign matrix dimensions to it by \code{dim(mat) <- c(nr, nc)}
    where \code{nr, nc} are the numbers of rows and columns. The
    resulting object \code{mat} is both a factor and a vector.
    \item
    Supply \code{mat} as a one-dimensional factor
    and specify the arguments \code{xcol} and \code{yrow}
    to determine the dimensions of the image.
    \item
    Use the functions
    \code{\link{cut.im}} or \code{\link{eval.im}} to make factor-valued
    images from other images).
  }
    
  For a description of the methods available for pixel image objects,
  see \code{\link{im.object}}.

  To convert other kinds of data to a pixel image (for example,
  functions or windows), use \code{\link{as.im}}.
}
\seealso{
  \code{\link{im.object}},
  \code{\link{as.im}},
  \code{\link{as.matrix.im}},
  \code{\link{[.im}},
  \code{\link{eval.im}}
}
\section{Warnings}{
  The internal representation of images is likely to change in future
  releases of \pkg{spatstat}. The safe way to extract pixel values
  from an image object is to use \code{\link{as.matrix.im}}
  or \code{\link{[.im}}.
}
\examples{
   vec <- rnorm(1200)
   mat <- matrix(vec, nrow=30, ncol=40)
   whitenoise <- im(mat)
   whitenoise <- im(mat, xrange=c(0,1), yrange=c(0,1))
   whitenoise <- im(mat, xcol=seq(0,1,length=40), yrow=seq(0,1,length=30))
   whitenoise <- im(vec, xcol=seq(0,1,length=40), yrow=seq(0,1,length=30))
   plot(whitenoise)

   # Factor-valued images:
   f <- factor(letters[1:12])
   dim(f) <- c(3,4)
   Z <- im(f)

   # Factor image from other image:
   cutwhite <- cut(whitenoise, 3)
   plot(cutwhite)

   # Factor image from raw data
   cutmat <- cut(mat, 3)
   dim(cutmat) <- c(30,40)
   cutwhite <- im(cutmat)
   plot(cutwhite)
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{manip}
\keyword{datagen}
 
 
