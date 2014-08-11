\name{im.object}
\alias{im.object}
\title{Class of Images}
\description{
  A class \code{"im"} to represent a two-dimensional pixel image.
}
\details{
  An object of this class represents
  a two-dimensional pixel image. It specifies
  \itemize{
    \item the dimensions of the rectangular array of pixels
    \item \eqn{x} and \eqn{y} coordinates for the pixels
    \item a numeric value (``grey value'') at each pixel
  }
  If \code{X} is an object of type \code{im},
  it contains the following elements:
  \tabular{ll}{
    \code{v} \tab matrix of values \cr
    \code{dim} \tab dimensions of matrix \code{v} \cr
    \code{xrange} \tab range of \eqn{x} coordinates of image window \cr
    \code{yrange} \tab range of \eqn{y} coordinates of image window \cr
    \code{xstep} \tab width of one pixel \cr
    \code{ystep} \tab height of one pixel \cr
    \code{xcol} \tab vector of \eqn{x} coordinates of centres of pixels \cr
    \code{yrow} \tab vector of \eqn{y} coordinates of centres of pixels 
  }
  Users are strongly advised not to manipulate these entries
  directly.

  Objects of class \code{"im"}
  may be created by the functions
  \code{\link{im}} and \code{\link{as.im}}.
  Image objects are also returned by various functions including
  \code{\link{distmap}}, \code{\link{Kmeasure}}, \code{\link{setcov}}, 
  \code{\link{eval.im}} and \code{\link{cut.im}}.

  Image objects may be displayed using the methods
  \code{\link{plot.im}}, \code{image.im}, \code{\link{persp.im}}
  and \code{contour.im}. There are also methods 
  \code{\link{print.im}} for printing information about an image,
  \code{\link{summary.im}} for summarising an image,
  \code{\link{mean.im}} for calculating the average pixel value, 
  \code{\link{hist.im}} for plotting a histogram of pixel values,
  \code{\link{quantile.im}} for calculating quantiles of pixel values,
  and \code{\link{cut.im}} for dividing the range of pixel values into
  categories.

  Pixel values in an image may be extracted
  using the subset operator \code{\link{[.im}}.
  To extract all pixel values from an image object,
  use \code{\link{as.matrix.im}}.
  The levels of a factor-valued image can be extracted and
  changed with \code{levels} and \code{levels<-}.

  Calculations involving one or more images (for example,
  squaring all the pixel values in an image, converting numbers to 
  factor levels, or 
  subtracting one image from another) can often be done
  easily using \code{\link{eval.im}}.
  To find all pixels satisfying
  a certain constraint, use \code{\link{solutionset}}.
  
  Note carefully that the entry \code{v[i,j]}
  gives the pixel value at the location \code{(xcol[j],yrow[i]}.
  That is, the \bold{row} index of the matrix \code{v} corresponds
  to increasing \bold{y} coordinate, while the column index of \code{mat}
  corresponds to increasing \bold{x} coordinate.
  Thus \code{yrow} has one entry for each row of \code{v}
  and \code{xcol} has one entry for each column of \code{v}.
  Under the usual convention in \R, a correct
  display of the image would be obtained by transposing the matrix, e.g.
  \code{image.default(xcol, yrow, t(v))}, if you wanted to do it by hand.
}
\seealso{
  \code{\link{im}},
  \code{\link{as.im}},
  \code{\link{plot.im}},
  \code{\link{persp.im}},
  \code{\link{eval.im}},
  \code{\link{[.im}}
}
\section{Warnings}{
  The internal representation of images is likely to change in future
  releases of \pkg{spatstat}. Do not address the entries in an
  image directly. To extract all pixel values from an image object,
  use \code{\link{as.matrix.im}}.
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{attribute}
 
 
