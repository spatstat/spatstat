\name{increment.fv}
\alias{increment.fv}
\title{
  Increments of a Function
}
\description{
  Compute the change in the value of a function \code{f}
  when the function argument increases by \code{delta}.
}
\usage{
increment.fv(f, delta)
}
\arguments{
  \item{f}{
    Object of class \code{"fv"} representing a function.
  }
  \item{delta}{
    Numeric. The increase in the value of the function argument.
  }
}
\details{
  This command computes the new function
  \deqn{g(x) = f(x+h) - f(x-h)}
  where \code{h = delta/2}. The value of \eqn{g(x)} is
  the change in the value of \eqn{f} over an interval of length
  \code{delta} centred at \eqn{x}.
}
\value{
  Another object of class \code{"fv"} compatible with \code{X}.
}
\author{Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  Rolf Turner \email{r.turner@auckland.ac.nz}
  
  and Ege Rubak \email{rubak@math.aau.dk}
  
}
\seealso{
  \code{\link{fv.object}}, 
  \code{\link{deriv.fv}}
}
\examples{
  plot(increment.fv(Kest(cells), 0.05))
}
\keyword{spatial}
\keyword{math}
\keyword{nonparametric}
