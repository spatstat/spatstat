\name{pool.rat}
\alias{pool.rat}
\title{
  Pool Data from Several Ratio Objects
}
\description{
  Pool the data from several ratio objects
  (objects of class \code{"rat"})
  and compute a pooled estimate.
}
\usage{
\method{pool}{rat}(...)
}
\arguments{
  \item{\dots}{
    Objects of class \code{"rat"}.
  }
}
\details{
  The function \code{\link{pool}} is generic. This is the method for the
  class \code{"rat"} of ratio objects. It is used to
  combine several estimates of the same quantity
  when each estimate is a ratio. 

  Each of the arguments \code{\dots} must be an object of class
  \code{"rat"} representing a ratio object (basically a
  numerator and a denominator; see \code{\link{rat}}).
  We assume that these ratios are all estimates of the same quantity.

  If the objects are called \eqn{R_1, \ldots, R_n}{R[1], \dots, R[n]}
  and if \eqn{R_i}{R[i]} has numerator \eqn{Y_i}{Y[i]} and
  denominator \eqn{X_i}{X[i]}, so that notionally
  \eqn{R_i = Y_i/X_i}{R[i] = Y[i]/X[i]}, then the pooled estimate is the
  ratio-of-sums estimator
  \deqn{
    R = \frac{\sum_i Y_i}{\sum_i X_i}.
  }{
    R = (Y[1]+\dots+Y[n])/(X[1]+\dots+X[n]).
  }
  The standard error of \eqn{R} is computed using the delta method
  as described in Baddeley \emph{et al.} (1993)
  or Cochran (1977, pp 154, 161).

  This calculation is implemented only for certain classes of objects
  where the arithmetic can be performed.
  
  This calculation is currently implemented only for objects which
  also belong to the class \code{"fv"} (function value tables).
  For example, if \code{\link{Kest}} is called with argument
  \code{ratio=TRUE}, the result is a suitable object (belonging to the classes
  \code{"rat"} and \code{"fv"}).

  Warnings or errors will be issued if the ratio objects \code{\dots}
  appear to be incompatible. However, the code is not smart enough to
  decide whether it is sensible to pool the data.
}
\value{
  An object of the same class as the input.
}
\seealso{
  \code{\link{rat}},
  \code{\link{pool}},
  \code{\link{Kest}}
}
\examples{
   K1 <- Kest(runifpoint(42), ratio=TRUE, correction="iso")   
   K2 <- Kest(runifpoint(42), ratio=TRUE, correction="iso")   
   K3 <- Kest(runifpoint(42), ratio=TRUE, correction="iso")
   K <- pool(K1, K2, K3)
   plot(K, pooliso ~ r, shade=c("hiiso", "loiso"))
}
\references{
  Baddeley, A.J, Moyeed, R.A., Howard, C.V. and Boyde, A. (1993)
  Analysis of a three-dimensional point pattern with replication.
  \emph{Applied Statistics} \bold{42}, 641--668.

  Cochran, W.G. (1977) 
  \emph{Sampling techniques}, 3rd edition.
  New York: John Wiley and Sons.
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{nonparametric}
