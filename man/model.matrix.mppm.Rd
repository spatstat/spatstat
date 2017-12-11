\name{model.matrix.mppm}
\alias{model.matrix.mppm}
\title{Extract Design Matrix of Point Process Model for Several Point Patterns}
\description{
  Given a point process model fitted to a list of point patterns,
  this function extracts the design matrix.
}
\usage{
   \method{model.matrix}{mppm}(object, ..., keepNA=TRUE, separate=FALSE)
}
\arguments{
  \item{object}{
    A point process model fitted to several point patterns. 
    An object of class \code{"mppm"}.
  }
  \item{\dots}{
    Other arguments (such as \code{na.action}) passed to
    \code{\link{model.matrix.lm}}.
  }
  \item{keepNA}{
    Logical. Determines whether rows containing \code{NA} values will be
    deleted or retained.
  }
  \item{separate}{
    Logical value indicating whether to split the model matrix into
    sub-matrices corresponding to each of the original point patterns.
  }
}
\details{
  This command is a method for the generic function
  \code{\link{model.matrix}}. It extracts the design matrix of a
  point process model fitted to several point patterns. 
  
  The argument \code{object} must be a fitted point process model
  (object of class \code{"mppm"}) produced by the 
  fitting algorithm \code{\link{mppm}}). This represents a
  point process model that has been fitted
  to a list of several point pattern datasets. See \code{\link{mppm}}
  for information.

  The result is a matrix with one column for every
  constructed covariate in the model, and one row
  for every quadrature point.
  
  If \code{separate=TRUE} this matrix will be split into
  sub-matrices corresponding to the original point patterns,
  and the result will be a list containing these matrices.
}
\value{
  A matrix (or list of matrices).
  Columns of the matrix are canonical covariates in the model.
}
\author{
  \spatstatAuthors.
}
\seealso{
  \code{\link{model.matrix}},
  \code{\link{mppm}}.
}
\examples{
   fit <- mppm(Points ~ Image + x, demohyper)
   head(model.matrix(fit))
   # matrix with three columns: '(Intercept)', 'x' and 'Image'
}
\keyword{spatial}
\keyword{models}
