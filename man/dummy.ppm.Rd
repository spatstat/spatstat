\name{dummy.ppm}
\alias{dummy.ppm}
\title{Extract Dummy Points Used to Fit a Point Process Model}
\description{
  Given a fitted point process model,
  this function extracts the `dummy points' of the
  quadrature scheme used to fit the model.
}
\usage{
  dummy.ppm(object, drop=FALSE)
}
\arguments{
  \item{object}{
    fitted point process model (an object of class \code{"ppm"}).
  }
  \item{drop}{
    Logical value determining whether to delete dummy points
    that were not used to fit the model.
  }
}
\value{
  A point pattern (object of class \code{"ppp"}).
}
\details{
  An object of class \code{"ppm"} represents a point process model
  that has been fitted to data. It is typically produced by
  the model-fitting algorithm \code{\link{ppm}}.
  
  The maximum pseudolikelihood algorithm in \code{\link{ppm}}
  approximates the pseudolikelihood
  integral by a sum over a finite set of quadrature points,
  which is constructed by augmenting
  the original data point pattern by a set of ``dummy'' points.
  The fitted model object returned by \code{\link{ppm}}
  contains complete information about this quadrature scheme.
  See \code{\link{ppm}} or \code{\link{ppm.object}} for further
  information.
  
  This function \code{dummy.ppm}
  extracts the dummy points of the quadrature scheme.
  A typical use of this function would be to count the number of dummy
  points, to gauge the accuracy of the approximation to the
  exact pseudolikelihood. 

  It may happen that some dummy points are not actually used in
  fitting the model (typically because the value of a covariate is \code{NA}
  at these points). The argument \code{drop} specifies whether these
  unused dummy points shall be deleted (\code{drop=TRUE}) or
  retained (\code{drop=FALSE}) in the return value.

  See \code{\link{ppm.object}} for a list of all operations that can be
  performed on objects of class \code{"ppm"}.
}
\seealso{
  \code{\link{ppm.object}},
  \code{\link{ppp.object}},
  \code{\link{ppm}}
}
\examples{
 data(cells)
 fit <- ppm(cells, ~1, Strauss(r=0.1))
 X <- dummy.ppm(fit)
 npoints(X)
 # this is the number of dummy points in the quadrature scheme
}
\author{Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  and Rolf Turner \email{r.turner@auckland.ac.nz}
  
}
\keyword{spatial}
\keyword{utilities}
\keyword{models}
