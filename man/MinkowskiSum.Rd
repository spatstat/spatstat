\name{MinkowskiSum}
\alias{MinkowskiSum}
\alias{\%+\%}  %DoNotExport 
%NAMESPACE export("%+%")
\title{Minkowski Sum of Windows}
\description{
  Compute the Minkowski sum of two spatial windows.
}
\usage{
MinkowskiSum(A, B)

A \%+\% B
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
  The operator \code{A \%+\% B} and function \code{MinkowskiSum(A,B)}
  are synonymous: they both compute the
  Minkowski sum of the windows \code{A} and \code{B}.

  The Minkowski sum 
  of two spatial regions \eqn{A} and \eqn{B}
  is another region, formed by taking all possible pairs of points,
  one in \eqn{A} and one in \eqn{B}, and adding them as vectors.
  The Minkowski Sum \eqn{A \oplus B}{A \%+\% B}
  is the set of all points \eqn{a+b} where \eqn{a} is in \eqn{A}
  and \eqn{b} is in \eqn{B}.
  A few common facts about the Minkowski sum are:
  \itemize{
    \item
    The sum is symmetric:
    \eqn{A \oplus B = B \oplus A}{A \%+\% B = B \%+\% A}.
    \item
    If \eqn{B} is a single point, then \eqn{A \oplus B}{A \%+\% B}
    is a shifted copy of \eqn{A}.
    \item
    If \eqn{A} is a square of side length \eqn{a},
    and \eqn{B} is a square of side length \eqn{b},
    with sides that are parallel to the coordinate axes,
    then \eqn{A \oplus B}{A \%+\% B} is a square of side length \eqn{a+b}.
    \item
    If \eqn{A} and \eqn{B} are discs of radius \eqn{r} and \eqn{s}
    respectively, then \eqn{A \oplus B}{A \%+\% B} is a disc
    of redius \eqn{r+s}.
    \item
    If \eqn{B} is a disc of radius \eqn{r} centred at the origin,
    then \eqn{A \oplus B}{A \%+\% B} is equivalent to the
    \emph{morphological dilation} of \eqn{A} by distance \eqn{r}.
    See \code{\link{dilation}}.
  }

  The algorithm currently computes the Minkowski sum as a polygonal
  window using the \pkg{polyclip} library. 
  It will be quite slow if applied to binary mask windows.
}
\seealso{
  \code{\link{dilation}}
}
\examples{
  opa <- par(mfrow=c(1,2))
  FR <- grow.rectangle(Frame(letterR), 0.3)
  plot(FR, main="")
  plot(letterR, add=TRUE, lwd=2, hatch=TRUE, hatchargs=list(texture=5))
  B <- square(0.2)
  plot(shift(B, vec=c(3.675, 3)),
       add=TRUE, border="red", lwd=2)
  plot(FR, main="")
  plot(letterR, add=TRUE)  
  plot(letterR \%+\% B, add=TRUE, border="blue", lwd=2,
         hatch=TRUE, hatchargs=list(col="blue"))
  par(opa)
}
\author{
  Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
}
\keyword{spatial}
\keyword{math}
