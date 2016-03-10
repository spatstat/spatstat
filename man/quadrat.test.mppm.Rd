\name{quadrat.test.mppm}
\alias{quadrat.test.mppm}
\title{Chi-Squared Test for Multiple Point Process Model Based on
  Quadrat Counts}
\description{
  Performs a chi-squared goodness-of-fit test of a 
  Poisson point process model fitted to multiple point patterns.
}
\usage{
  \method{quadrat.test}{mppm}(X, ...)
}
\arguments{
  \item{X}{
    An object of class \code{"mppm"} representing a
    point process model fitted to multiple point patterns.
    It should be a Poisson model.
  }
  \item{\dots}{
    Arguments passed to \code{\link{quadrat.test.ppm}}
    which determine the size of the quadrats.
  }
}
\details{
  This function performs a \eqn{\chi^2}{chi^2} test of goodness-of-fit
  for a Poisson point process model, based on quadrat counts.
  It can also be used to perform a test of Complete Spatial Randomness
  for a list of point patterns.

  The function \code{quadrat.test} is generic, with methods for
  point patterns (class \code{"ppp"}), point process models
  (class \code{"ppm"}) and
  multiple point process models (class
  \code{"mppm"}).

  For this function, the argument \code{X} should be a
  multiple point process model (object of class \code{"mppm"})
  obtained by fitting a point process model to a list of
  point patterns using the function \code{\link{mppm}}.

  To perform the test, the data point patterns are extracted from \code{X}.
  For each point pattern
  \itemize{
    \item the window of observation is divided
    into rectangular tiles, and the number of data points in each tile is
    counted, as described in \code{\link{quadratcount}}.
    \item 
    The expected number of points in each quadrat is calculated,
    as determined by the fitted model.
  }
  Then we perform a single \eqn{\chi^2}{chi^2} test of goodness-of-fit
  based on these observed and expected counts. 
}
\section{Testing Complete Spatial Randomness}{
  If the intention is to test Complete Spatial Randomness (CSR)
  there are two options:
  \itemize{
    \item CSR with the same intensity of points in each point pattern;
    \item CSR with a different, unrelated intensity of points in each
    point pattern.
  }
  In the first case, 
  suppose \code{P} is a list of point patterns we want to test.
  Then fit the multiple model \code{fit1 <- mppm(P, ~1)} which signifies a
  Poisson point process model with a constant intensity. Then
  apply \code{quadrat.test(fit1)}.

  In the second case, fit the model code{fit2 <- mppm(P, ~id)}
  which signifies a Poisson point process with a different constant
  intensity for each point pattern. Then apply \code{quadrat.test(fit2)}.
}

\value{
  An object of class \code{"htest"}.
  Printing the object gives comprehensible output
  about the outcome of the test.
  The \eqn{p}-value of the test is stored in the
  component \code{p.value}.

  The return value also belongs to
  the special class \code{"quadrat.test"}. Plotting the object
  will display, for each window, the position of the quadrats,
  annotated by their observed and expected
  counts and the Pearson residuals. See the examples.

  The return value also has an attribute \code{"components"}
  which is a list containing the results of 
  \eqn{\chi^2}{chi^2} tests of goodness-of-fit
  for each individual point pattern.
}
\seealso{
  \code{\link{mppm}},
  \code{\link{quadrat.test}}
}
\examples{
  H <- hyperframe(X=waterstriders)
  # Poisson with constant intensity for all patterns
  fit1 <- mppm(X~1, H)
  quadrat.test(fit1, nx=2)

  # uniform Poisson with different intensity for each pattern
  fit2 <- mppm(X ~ id, H)
  quadrat.test(fit2, nx=2)
}
\references{
  Baddeley, A., Rubak, E. and Turner, R. (2015)
  \emph{Spatial Point Patterns: Methodology and Applications with R}.
  London: Chapman and Hall/CRC Press. 
}
\author{
  Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}, Ida-Maria Sintorn and Leanne Bischoff.
  Implemented by 
  Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  Rolf Turner \email{r.turner@auckland.ac.nz}
  
  and Ege Rubak \email{rubak@math.aau.dk}
  
}
\keyword{spatial}
\keyword{htest}

