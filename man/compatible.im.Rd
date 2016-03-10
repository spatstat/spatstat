\name{compatible.im}  %DontDeclareMethods
\alias{compatible.im}
\title{Test Whether Pixel Images Are Compatible}
\description{
  Tests whether two or more pixel image objects have compatible dimensions.
}
\usage{
 \method{compatible}{im}(A, B, \dots, tol=1e-6)
}
\arguments{
  \item{A,B,\dots}{Two or more pixel images (objects of class \code{"im"}).}
  \item{tol}{Tolerance factor}
}
\details{
  This function tests whether the pixel images \code{A} and \code{B}
  (and any additional images \code{\dots})
  have compatible pixel dimensions. They are compatible if they have
  the same number of rows and columns, the same physical pixel
  dimensions, and occupy the same rectangle in the plane.

  The argument \code{tol} specifies the maximum tolerated error
  in the pixel coordinates, expressed as a
  fraction of the dimensions of a single pixel.
}
\value{
  Logical value: \code{TRUE} if the images are compatible, and \code{FALSE}
  if they are not.
}
\seealso{
  \code{\link{eval.im}},
  \code{\link{harmonise.im}},
  \code{\link{commonGrid}}
}
\author{Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  and Rolf Turner \email{r.turner@auckland.ac.nz}
  
}
\keyword{spatial}
\keyword{manip}
