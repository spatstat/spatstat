\name{dffit.ppm}
\alias{dffit}
\alias{dffit.ppm}
\title{
  Case Deletion Effect Measure of Fitted Model
}
\description{
  Computes the case deletion effect measure \code{DFFIT}
  for a fitted model.
}
\usage{
dffit(object, \dots)

\method{dffit}{ppm}(object, \dots, collapse = FALSE, dfb = NULL)
}
\arguments{
  \item{object}{
    A fitted model, such as a point process model
    (object of class \code{"ppm"}).
  }
  \item{\dots}{
    Additional arguments passed to \code{\link{dfbetas.ppm}}.
  }
  \item{collapse}{
    Logical value specifying whether to collapse the vector-valued
    measure to a scalar-valued measure by adding all the components.
  }
  \item{dfb}{
    Optional. The result of \code{dfbetas(object)},
    if it has already been computed.
  }
}
\details{
  The case deletion effect measure \code{DFFIT} is a model diagnostic
  traditionally used for regression models. In that context,
  \code{DFFIT[i,j]} is the negative change, in the value of the
  \code{j}th term in the linear predictor, that would occur if the \code{i}th
  data value was deleted. It is closely related to the
  diagnostic \code{DFBETA}.

  For a spatial point process model, \code{dffit} computes
  the analogous spatial case deletion diagnostic, described in
  Baddeley, Rubak and Turner (2018).
}
\value{
  A measure (object of class \code{"msr"}).
}
\references{
  Baddeley, A., Rubak, E. and Turner, R. (2018)
  \emph{Leverage and influence diagnostics for Gibbs spatial point
    processes}.
  In preparation.
}
\author{
  \spatstatAuthors.
}
\seealso{
  \code{\link{dfbetas.ppm}}
}
\examples{
   \testonly{op <- spatstat.options(npixel=32)}
   X <- rpoispp(function(x,y) { exp(3+3*x) })
   fit <- ppm(X ~x+y)
   \testonly{fit <- ppm(X ~x+y, nd=16)}
   plot(dffit(fit))
   plot(dffit(fit, collapse=TRUE))
   \testonly{spatstat.options(op)}
   
}
\keyword{spatial}
\keyword{models}
