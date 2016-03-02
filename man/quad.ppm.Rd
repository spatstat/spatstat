\name{quad.ppm}
\alias{quad.ppm}
\title{Extract Quadrature Scheme Used to Fit a Point Process Model}
\description{
  Given a fitted point process model,
  this function extracts the 
  quadrature scheme used to fit the model.
}
\usage{
  quad.ppm(object, drop=FALSE, clip=FALSE)
}
\arguments{
  \item{object}{
    fitted point process model (an object of class \code{"ppm"}
    or \code{"kppm"} or \code{"lppm"}).
  }
  \item{drop}{
    Logical value determining whether to delete quadrature points
    that were not used to fit the model.
  }
  \item{clip}{
    Logical value determining whether to erode the window,
    if \code{object} was fitted using the border correction.
    See Details.
  }
}
\value{
  A quadrature scheme (object of class \code{"quad"}).
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
  
  This function \code{quad.ppm} extracts the quadrature scheme.
  A typical use of this function would be to inspect the quadrature scheme
  (points and weights) to gauge the accuracy of the approximation to the
  exact pseudolikelihood.

  Some quadrature points may not have been used in
  fitting the model. This happens if the border correction is used,
  and in other cases (e.g. when the value of a covariate is \code{NA}
  at these points). The argument \code{drop} specifies whether these
  unused quadrature points shall be deleted (\code{drop=TRUE}) or
  retained (\code{drop=FALSE}) in the return value.

  The quadrature scheme has a \emph{window}, which by default is set to
  equal the window of the original data. However this window may be
  larger than the actual domain of integration of the pseudolikelihood
  or composite likelihood that was used to fit the model.
  If \code{clip=TRUE} then the window of the quadrature scheme is
  set to the actual domain of integration. This option only has an effect
  when the model was fitted using the border correction; then
  the window is obtained by eroding the original data window
  by the border correction distance. 
  
  See \code{\link{ppm.object}} for a list of all operations that can be
  performed on objects of class \code{"ppm"}.
  See \code{\link{quad.object}} for a list of all operations that can be
  performed on objects of class \code{"quad"}.

  This function can also be applied to objects of class \code{"kppm"}
  and \code{"lppm"}.
}
\seealso{
  \code{\link{ppm.object}},
  \code{\link{quad.object}},
  \code{\link{ppm}}
}
\examples{
 fit <- ppm(cells ~1, Strauss(r=0.1))
 Q <- quad.ppm(fit)
 \dontrun{plot(Q)}
 npoints(Q$data)
 npoints(Q$dummy)
}
\author{Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  and Rolf Turner \email{r.turner@auckland.ac.nz}
  
}
\keyword{spatial}
\keyword{manip}
\keyword{models}
