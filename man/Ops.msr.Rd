\name{Ops.msr}
\alias{Ops.msr}
\title{Arithmetic Operations on Measures}
\description{
  These group generic methods for the class \code{"msr"}
  allow the arithmetic operators
  \code{+}, \code{-}, \code{*} and \code{/} to be applied
  directly to measures. 
}
\usage{
## S3 methods for group generics have prototypes:
\special{Ops(e1, e2)}
%NAMESPACE S3method("Ops", "msr")
}
\arguments{
  \item{e1, e2}{objects of class \code{"msr"}.}
}
\details{
  Arithmetic operators on a measure \code{A} are only defined
  in some cases. The arithmetic operator is effectively applied
  to the value of \code{A(W)} for every spatial domain \code{W}.
  If the result is a measure, then this operation is valid.
  
  If \code{A} is a measure (object of class \code{"msr"})
  then the operations \code{-A} and \code{+A} are defined.
  
  If \code{A} and \code{B} are measures,
  then \code{A + B} and \code{A - B} are defined.
  (This is implemented only when \code{A} and \code{B}
  are represented using the same quadrature scheme.)

  If \code{A} is a measure and \code{z} is a numeric value,
  then \code{A * z} and \code{A / z} are defined,
  and \code{z * A} is defined.
}
\seealso{
  \code{\link{with.msr}}
}
\examples{
   X <- rpoispp(function(x,y) { exp(3+3*x) })
   fit <- ppm(X, ~x+y)
   rp <- residuals(fit, type="pearson")
   rp

   -rp
   2 * rp
   rp /2

   rp - rp

   rr <- residuals(fit, type="raw")
   rp - rr
}
\author{
  \spatstatAuthors.
}
\keyword{spatial}
\keyword{methods}
