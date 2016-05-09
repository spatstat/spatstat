\name{erosionAny}
\alias{erosionAny}
\alias{\%(-)\%}  %DoNotExport 
%NAMESPACE export("%(-)%")
\title{Morphological Erosion of Windows}
\description{
  Compute the morphological erosion of one spatial window by another.
}
\usage{
erosionAny(A, B)

A \%(-)\% B
}
\arguments{
  \item{A,B}{
    Windows (objects of class \code{"owin"}).
  }
}
\value{
  Another window (object of class \code{"owin"}).
}
\details{
  The operator \code{A \%(-)\% B} and function \code{erosionAny(A,B)}
  are synonymous: they both compute the
  morphological erosion of the window \code{A} by the window \code{B}.

  The morphological erosion
  \eqn{A \ominus B}{A \%(-)\% B}
  of region \eqn{A} by region \eqn{B}
  is the spatial region consisting of all vectors \eqn{z}
  such that, when \eqn{B} is shifted by the vector \eqn{z}, the result
  is a subset of \eqn{A}.

  Equivalently
  \deqn{
    A \ominus B = ((A^c \oplus (-B))^c
  }{
    (A^c \%+\% (-B))^c
  }
  where \eqn{\oplus}{\%+\%} is the Minkowski sum,
  \eqn{A^c} denotes the set complement, and \eqn{(-B)} is
  the reflection of \eqn{B} through the origin, consisting of all
  vectors \eqn{-b} where \eqn{b} is a point in \eqn{B}.

  If \code{B} is a disc of radius \code{r}, then
  \code{erosionAny(A, B)} is equivalent to \code{erosion(A, r)}.
  See \code{\link{erosion}}.
  
  The algorithm currently computes the result as a polygonal
  window using the \pkg{polyclip} library. 
  It will be quite slow if applied to binary mask windows.
}
\seealso{
  \code{\link{erosion}},
  \code{\link{MinkowskiSum}}
}
\examples{
  B <- square(c(-0.1, 0.1))
  RminusB <- letterR \%(-)\% B
  FR <- grow.rectangle(Frame(letterR), 0.3)
  plot(FR, main="", type="n")
  plot(letterR, add=TRUE, lwd=2, hatch=TRUE, box=FALSE)
  plot(RminusB, add=TRUE, col="blue", box=FALSE)
  plot(shift(B, vec=c(3.49, 2.98)),
       add=TRUE, border="red", lwd=2)
}
\author{
  Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au},
 Rolf Turner \email{r.turner@auckland.ac.nz}
 and Ege Rubak \email{rubak@math.aau.dk}
}
\keyword{spatial}
\keyword{math}
