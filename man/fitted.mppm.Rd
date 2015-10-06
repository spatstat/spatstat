\name{fitted.mppm}
\alias{fitted.mppm}
\title{Fitted Conditional Intensity for Multiple Point Process Model}
\description{
  Given a point process model fitted to multiple point patterns, compute the
  fitted conditional intensity of the model at the points of each data
  pattern, or at the points of the quadrature schemes used to fit the
  model.
}
\usage{
 \method{fitted}{mppm}(object, ..., type = "lambda", dataonly = FALSE)
}
\arguments{
  \item{object}{
    The fitted model. An object of class \code{"mppm"}
    obtained from \code{\link{mppm}}.
  }
  \item{\dots}{Ignored.}
  \item{type}{
    Type of fitted values: either \code{"trend"} for the spatial trend,
    or \code{"lambda"} or \code{"cif"} for the conditional intensity.
  }
  \item{dataonly}{
    If \code{TRUE}, fitted values are computed only for the points
    of the data point patterns. If \code{FALSE}, fitted values are
    computed for the points of the quadrature schemes used to fit the
    model.
  }
}
\details{
  This function evaluates the conditional intensity
  \eqn{\hat\lambda(u,x)}{lambdahat(u,x)}
  or spatial trend \eqn{\hat{b(u)}}{bhat(u)}
  of the fitted point process model for
  certain locations \eqn{u}, for each of the original point patterns \eqn{x}
  to which the model was fitted.

  The locations \eqn{u} at which the fitted conditional intensity/trend is
  evaluated, are the points of the quadrature schemes used to fit the
  model in \code{\link{mppm}}. They include the data points (the points of the
  original point pattern datasets) and other ``dummy'' points  in
  the window of observation.

  Use \code{\link{predict.mppm}} to compute the fitted conditional intensity at
  other locations or with other values of the explanatory variables.
}
\value{
  A list of vectors (one for each row of the original hyperframe,
  i.e. one vector for each of the original point patterns)
  containing the values of the fitted conditional intensity
  or (if \code{type="trend"}) the fitted spatial trend.

  Entries in these vector correspond to the quadrature points (data
  or dummy points) used to fit the model. The quadrature points can
  be extracted from \code{object} by \code{\link{quad.mppm}(object)}.
}
\references{
  Baddeley, A., Rubak, E. and Turner, R. (2015)
  \emph{Spatial Point Patterns: Methodology and Applications with R}.
  London: Chapman and Hall/CRC Press. 
}
\author{
  Adrian Baddeley, Ida-Maria Sintorn and Leanne Bischoff.
  Implemented by 
  Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  Rolf Turner
  \email{r.turner@auckland.ac.nz}
  and Ege Rubak
  \email{rubak@math.aau.dk}
}
\examples{
   model <- mppm(Bugs ~ x, data=hyperframe(Bugs=waterstriders),
                           interaction=Strauss(7))
   cifs <- fitted(model)
}
\seealso{
  \code{\link{mppm}},
  \code{\link{predict.mppm}}
}
\keyword{spatial}
\keyword{models}
