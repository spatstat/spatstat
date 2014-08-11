\name{model.matrix.ppm}
\alias{model.matrix.ppm}
\alias{model.matrix.kppm}
\alias{model.matrix.lppm}
\title{Extract Design Matrix from Point Process Model}
\description{
  Given a point process model that has been
  fitted to spatial point pattern data,
  this function extracts the design matrix of the model.
}
\usage{
   \method{model.matrix}{ppm}(object, data=model.frame(object), ..., keepNA=TRUE)

   \method{model.matrix}{kppm}(object, data=model.frame(object), ...,
   keepNA=TRUE)

   \method{model.matrix}{lppm}(object, data=model.frame(object), ..., keepNA=TRUE)
}
\arguments{
  \item{object}{
    The fitted point process model. An object of class \code{"ppm"}
    or \code{"kppm"} or \code{"lppm"}.
  }
  \item{data}{
    A model frame, containing the data required for the Berman-Turner device.
  }
  \item{keepNA}{
    Logical. Determines whether rows containing NA values will be
    deleted or retained.
  }
  \item{\dots}{
    Other arguments (such as \code{na.action}) passed to
    \code{\link{model.matrix.lm}}.
  }
}
\details{
  These commands are methods for the generic function
  \code{\link{model.matrix}}. They extracts the design matrix of a
  spatial point process model (class \code{"ppm"} or \code{"kppm"}
  or \code{"lppm"}).

  More precisely, this command extracts
  the design matrix of the generalised linear model associated with
  a spatial point process model. 
  
  The \code{object} must be a fitted point process model
  (object of class \code{"ppm"} or \code{"kppm"} or \code{"lppm"})
  fitted to spatial point pattern data.
  Such objects are produced by the model-fitting
  functions \code{\link{ppm}}, \code{\link{kppm}} and \code{\link{lppm}}.

  The methods \code{model.matrix.ppm}, \code{model.matrix.kppm}
  and \code{model.matrix.lppm} 
  extract the model matrix for the GLM.

  The result is a matrix, with one row for every quadrature point
  in the fitting procedure, and one column for every constructed
  covariate in the design matrix.

  If there are \code{NA} values in the covariates,
  the argument \code{keepNA} determines whether to retain or delete
  the corresponding rows of the model matrix. The default
  \code{keepNA=TRUE} is to retain them. Note that this differs from
  the default behaviour of many other methods for \code{model.matrix},
  which typically delete rows containing \code{NA}.
  
  The quadrature points themselves can be extracted using
  \code{\link{quad.ppm}}.
}
\value{
  A matrix. Columns of the matrix are covariates in the model.
  Rows of the matrix correspond to quadrature points
  in the fitting procedure (provided \code{keepNA=TRUE}).
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\seealso{
  \code{\link{model.matrix}},
  \code{\link{model.images}},
  \code{\link{ppm}},
  \code{\link{kppm}},
  \code{\link{lppm}},
  \code{\link{ppm.object}},
  \code{\link{quad.ppm}},
  \code{\link{residuals.ppm}}
}

\examples{
   fit <- ppm(cells, ~x)
   head(model.matrix(fit))
   # matrix with two columns: '(Intercept)' and 'x'
   kfit <- kppm(redwood, ~x, "Thomas")
   m <- model.matrix(kfit)
}
\keyword{spatial}
\keyword{models}
