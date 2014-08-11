\name{eem}
\alias{eem}
\title{
  Exponential Energy Marks
}
\description{
  Given a point process model fitted to a point pattern,
  compute the Stoyan-Grabarnik diagnostic ``exponential energy marks''
  for the data points. 
}
\usage{
  eem(fit, check=TRUE)
}
\arguments{
  \item{fit}{
    The fitted point process model. An object of class \code{"ppm"}.
  }
  \item{check}{
    Logical value indicating whether to check the internal format
    of \code{fit}. If there is any possibility that this object
    has been restored from a dump file, or has otherwise lost track of
    the environment where it was originally computed, set
    \code{check=TRUE}. 
  }
}
\value{
  A vector containing the values of the exponential energy mark
  for each point in the pattern.
}
\details{
  Stoyan and Grabarnik (1991) proposed a diagnostic
  tool for point process models fitted to spatial point pattern data.
  Each point \eqn{x[i]}{x_i} of the data pattern \eqn{X}
  is given a `mark' or `weight'
  \deqn{m[i] = 1/lambda-hat(x[i],X)}{m_i = \frac 1 {\hat\lambda(x_i,X)}}
  where \eqn{lambda-hat(x[i],X)}{\hat\lambda(x_i,X)}
  is the conditional intensity of the fitted model.
  If the fitted model is correct, then the sum of these marks
  for all points in a region \eqn{B} has expected value equal to the
  area of \eqn{B}.
  
  The argument \code{fit} must be a fitted point process model
  (object of class \code{"ppm"}). Such objects are produced by the maximum
  pseudolikelihood fitting algorithm \code{\link{ppm}}).
  This fitted model object contains complete
  information about the original data pattern and the model that was
  fitted to it.

  The value returned by \code{eem} is the vector
  of weights \eqn{m[i]}{m_i} associated with the points \eqn{x[i]}{x_i}
  of the original data pattern. The original data pattern
  (in corresponding order) can be
  extracted from \code{fit} using \code{\link{data.ppm}}.
  
  The function \code{\link{diagnose.ppm}}
  produces a set of sensible diagnostic plots based on these weights.
}
\references{
  Stoyan, D. and Grabarnik, P. (1991)
  Second-order characteristics for stochastic structures connected with
  Gibbs point processes.
  \emph{Mathematische Nachrichten}, 151:95--100.
}
\seealso{
 \code{\link{diagnose.ppm}},
 \code{\link{ppm.object}},
 \code{\link{data.ppm}},
 \code{\link{residuals.ppm}},
 \code{\link{ppm}}
}
\examples{
    data(cells)
    fit <- ppm(cells, ~x, Strauss(r=0.15))
    ee <- eem(fit)
    sum(ee)/area.owin(cells$window) # should be about 1 if model is correct
    Y <- setmarks(cells, ee)
    plot(Y, main="Cells data\n Exponential energy marks")
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{models}
