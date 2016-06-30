\name{model.frame.ppm}
\alias{model.frame.ppm}
\alias{model.frame.kppm}
\alias{model.frame.dppm}
\alias{model.frame.lppm}
\title{
Extract the Variables in a Point Process Model
}
\description{
  Given a fitted point process model, this function
  returns a data frame containing all the variables needed to
  fit the model using the Berman-Turner device.
}
\usage{
 \method{model.frame}{ppm}(formula, ...)

 \method{model.frame}{kppm}(formula, ...)

 \method{model.frame}{dppm}(formula, ...)

 \method{model.frame}{lppm}(formula, ...)
}
\arguments{
  \item{formula}{
    A fitted point process model.
    An object of class \code{"ppm"} or \code{"kppm"} or \code{"dppm"} or \code{"lppm"}.
  }
  \item{\dots}{
    Additional arguments passed to \code{\link{model.frame.glm}}.
  }
}
\details{
  The function \code{\link{model.frame}} is generic.
  These functions are method for \code{\link{model.frame}}
  for fitted point process models (objects of class \code{"ppm"}
  or \code{"kppm"} or \code{"dppm"} or \code{"lppm"}).

  The first argument should be a fitted point process model;
  it has to be named \code{formula} for consistency with the generic
  function.

  The result is a data frame containing all the variables used in
  fitting the model. The data frame has one row for each quadrature point
  used in fitting the model. The quadrature scheme can be extracted using
  \code{\link{quad.ppm}}.  
}
\value{
  A \code{data.frame} containing all the variables used in the
  fitted model, plus additional variables specified in \code{\dots}.
  It has an additional attribute \code{"terms"} containing information
  about the model formula. For details see \code{\link{model.frame.glm}}.
}
\references{
  Baddeley, A. and Turner, R. (2000)
  Practical maximum pseudolikelihood for spatial point patterns.
  \emph{Australian and New Zealand Journal of Statistics}
  \bold{42}, 283--322.
}
\seealso{
  \code{\link{ppm}},
  \code{\link{kppm}},
  \code{\link{dppm}},
  \code{\link{lppm}},
  \code{\link{model.frame}},
  \code{\link{model.matrix.ppm}}
}
\examples{
  fit <- ppm(cells ~ x)
  mf <- model.frame(fit)
  kfit <- kppm(redwood ~ x, "Thomas")
  kmf <- model.frame(kfit)
}
\author{\adrian
  
  
  \rolf
  
  and \ege
  
}
\keyword{spatial}
\keyword{models}

