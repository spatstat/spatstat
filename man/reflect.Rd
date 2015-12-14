\name{reflect} %DontDeclareMethods
\alias{reflect}
\alias{reflect.im}
\alias{reflect.default}
\title{Reflect In Origin}
\description{
  Reflects a geometrical object through the origin.
}
\usage{
  reflect(X)

  \method{reflect}{im}(X)

  \method{reflect}{default}(X)
}
\arguments{
  \item{X}{Any suitable dataset representing a two-dimensional
    object, such as a point pattern (object of class \code{"ppp"}),
    or a window (object of class \code{"owin"}).}
}
\value{
  Another object of the same type, representing the
  result of reflection.
}
\details{
  The object \code{X} is reflected through the origin.
  That is, each point in \code{X} with coordinates
  \eqn{(x,y)} is mapped to the position \eqn{(-x, -y)}.

  This is equivalent to applying the affine transformation with matrix
  \code{diag(c(-1,-1))}. It is also equivalent to rotation about the origin
  by 180 degrees.

  The command \code{reflect} is generic, with a method for
  pixel images and a default method.
}
\seealso{
  \code{\link{affine}},
  \code{\link{flipxy}}
}
\examples{
  plot(reflect(as.im(letterR)))
  plot(reflect(letterR), add=TRUE)
}
\author{Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  and Rolf Turner \email{r.turner@auckland.ac.nz}
  
}
\keyword{spatial}
\keyword{math}
