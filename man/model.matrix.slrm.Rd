\name{model.matrix.slrm}
\alias{model.matrix.slrm}
\title{Extract Design Matrix from Spatial Logistic Regression Model}
\description{
  This function extracts the design matrix of a
  spatial logistic regression model.
}
\usage{
   \method{model.matrix}{slrm}(object, ..., keepNA=TRUE)
}
\arguments{
  \item{object}{
    A fitted spatial logistic regression model. An object of class 
    \code{"slrm"}.
  }
  \item{\dots}{
    Other arguments (such as \code{na.action}) passed to
    \code{\link{model.matrix.lm}}.
  }
  \item{keepNA}{
    Logical. Determines whether rows containing \code{NA} values will be
    deleted or retained.
  }
}
\details{
  This command is a method for the generic function
  \code{\link{model.matrix}}. It extracts the design matrix of a
  spatial logistic regression.
  
  The \code{object} must be a fitted
  spatial logistic regression 
  (object of class \code{"slrm"}).
  Such objects are produced by the model-fitting
  function \code{\link{slrm}}.

  Usually the result is a matrix with one column for every
  constructed covariate in the model, and one row
  for every pixel in the grid used to fit the model.

  If \code{object} was fitted using split pixels (by calling
  \code{\link{slrm}} using the argument \code{splitby}) then the
  matrix has one row for every pixel or half-pixel.
}
\value{
  A matrix. Columns of the matrix are canonical covariates in the model.
}
\author{\adrian
  
  
  and \rolf
  
}
\seealso{
  \code{\link{model.matrix}},
  \code{\link{model.images}},
  \code{\link{slrm}}.
}
\examples{
   fit <- slrm(japanesepines ~x)
   head(model.matrix(fit))
   # matrix with two columns: '(Intercept)' and 'x'
}
\keyword{spatial}
\keyword{models}
