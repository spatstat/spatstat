\name{transmat}
\alias{transmat}
\title{
  Convert Pixel Array Between Different Conventions
}
\description{
  This function provides a simple way to convert arrays of pixel data
  between different display conventions.
}
\usage{
transmat(m, from, to)
}
\arguments{
  \item{m}{
    A matrix.
  }
  \item{from,to}{
    Specifications of the spatial arrangement of the pixels.
    See Details.
  }
}
\details{
  Pixel images are handled by many different software packages.
  In virtually all of these, the pixel values are stored in a matrix,
  and are accessed using the row and column indices of the matrix.
  However, different pieces of software use different conventions for
  mapping the matrix indices \eqn{[i,j]} to the spatial coordinates
  \eqn{(x,y)}.

  \itemize{
    \item
    In the \emph{Cartesian} convention, the first matrix index \eqn{i}
    is associated with the first Cartesian coordinate \eqn{x}, 
    and \eqn{j} is associated with \eqn{y}. This convention is used in
    \code{\link[graphics]{image.default}}.
    \item
    In the \emph{European reading order} convention, a matrix is displayed
    in the spatial coordinate system as it would be printed in a page of text:
    \eqn{i} is effectively associated with the negative \eqn{y} coordinate, 
    and \eqn{j} is associated with \eqn{x}. This convention is used in some
    image file formats.
    \item
    In the \code{spatstat} convention, \eqn{i}
    is associated with the increasing \eqn{y} coordinate, 
    and \eqn{j} is associated with \eqn{x}. This is also used in some
    image file formats.
  }
  To convert between these conventions, use the function
  \code{transmat}. If a matrix \code{m} contains pixel image data
  that is correctly displayed by software that uses the Cartesian convention,
  and we wish to convert it to the European reading convention, we can type
  \code{mm <- transmat(m, from="Cartesian", to="European")}.
  The transformed matrix \code{mm} will then be correctly displayed by
  software that uses the European convention. 

  Each of the arguments \code{from} and \code{to} can be one of the names
  \code{"Cartesian"}, \code{"European"} or \code{"spatstat"} (partially matched)
  or it can be a list specifying another convention. For example
  \code{to=list(x="-i", y="-j")!} specifies that rows of the output matrix 
  are expected to be displayed as vertical columns in the plot, 
  starting at the right side of the plot, as in the traditional 
  Chinese, Japanese and Korean writing order.
}
\value{
  Another matrix obtained by rearranging the entries of \code{m}.
}
\author{
Adrian Baddeley
\email{Adrian.Baddeley@uwa.edu.au}
\url{http://www.maths.uwa.edu.au/~adrian/}
Rolf Turner
\email{r.turner@auckland.ac.nz}
and Ege Rubak
\email{rubak@math.aau.dk}
}
\examples{
  opa <- par(mfrow=c(1,2))
  # image in spatstat format
  Z <- bei.extra$elev
  plot(Z, main="plot.im", ribbon=FALSE)
  m <- as.matrix(Z)
  # convert matrix to format suitable for display by image.default
  Y <- transmat(m, from="spatstat", to="Cartesian")
  image(Y, asp=0.5, main="image.default", axes=FALSE)
  par(opa)
}
\keyword{spatial}
\keyword{hplot}
\keyword{manip}
