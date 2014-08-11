\name{bermantest}
\alias{bermantest}
\alias{bermantest.ppm}
\alias{bermantest.ppp}
\alias{bermantest.lppm}
\alias{bermantest.lpp}
\title{Berman's Tests for Point Process Model} 
\description{
  Tests the goodness-of-fit of a Poisson point process model
  using methods of Berman (1986).
}
\usage{
bermantest(...)

\method{bermantest}{ppp}(X, covariate,
                         which = c("Z1", "Z2"),
        alternative = c("two.sided", "less", "greater"), ...)

\method{bermantest}{ppm}(model, covariate,
                         which = c("Z1", "Z2"),
               alternative = c("two.sided", "less", "greater"), ...)

\method{bermantest}{lpp}(X, covariate,
                         which = c("Z1", "Z2"),
        alternative = c("two.sided", "less", "greater"), ...)

\method{bermantest}{lppm}(model, covariate,
                         which = c("Z1", "Z2"),
               alternative = c("two.sided", "less", "greater"), ...)
}
\arguments{
  \item{X}{
    A point pattern (object of class \code{"ppp"} or \code{"lpp"}).
  }
  \item{model}{
    A fitted point process model (object of class \code{"ppm"} or
    \code{"lppm"}).
  }
  \item{covariate}{
    The spatial covariate on which the test will be based.
    An image (object of class \code{"im"})
    or a function.
  }
  \item{which}{
    Character string specifying the choice of test.
  }
  \item{alternative}{
    Character string specifying the alternative hypothesis.
  }
  \item{\dots}{
    Ignored.
  }
}
\details{
  These functions perform a goodness-of-fit test of a Poisson point
  process model fitted to point pattern data. The observed distribution
  of the values of a spatial covariate at the data points,
  and the predicted distribution of the same values under the model,
  are compared using either of two test statistics
  \eqn{Z_1}{Z[1]} and \eqn{Z_2}{Z[2]} proposed by Berman (1986).

  The function \code{bermantest} is generic, with methods for
  point patterns (\code{"ppp"} or \code{"lpp"})
  and point process models (\code{"ppm"} or \code{"lppm"}).
  \itemize{
    \item 
    If \code{X} is a point pattern dataset (object of class
    \code{"ppp"} or \code{"lpp"}), then
    \code{bermantest(X, ...)} performs a goodness-of-fit test of the
    uniform Poisson point process (Complete Spatial Randomness, CSR)
    for this dataset.
    \item
    If \code{model} is a fitted point process model
    (object of class \code{"ppm"} or \code{"lppm"})
    then \code{bermantest(model, ...)} performs
    a test of goodness-of-fit for this fitted model. In this case,
    \code{model} should be a Poisson point process.
  }
  
  The test is performed by comparing the observed distribution
  of the values of a spatial covariate at the data points,
  and the predicted distribution of the same covariate under the model.
  Thus, you must nominate a spatial covariate for this test.
  
  The argument \code{covariate} should be either a \code{function(x,y)}
  or a pixel image (object of class \code{"im"} containing the values
  of a spatial function.
  If \code{covariate} is an image, it should have numeric values,
  and its domain should cover the observation window of the
  \code{model}. If \code{covariate} is a function, it should expect
  two arguments \code{x} and \code{y} which are vectors of coordinates,
  and it should return a numeric vector of the same length
  as \code{x} and \code{y}.  

  First the original data point pattern is extracted from \code{model}.
  The values of the \code{covariate} at these data points are
  collected.

  Next the values of the \code{covariate} at all locations in the
  observation window are evaluated. The point process intensity
  of the fitted model is also evaluated at all locations in the window.

  \itemize{
    \item If \code{which="Z1"},
    the test statistic \eqn{Z_1}{Z[1]} is computed as follows.
    The sum \eqn{S} of the covariate values at all data
    points is evaluated. The predicted mean \eqn{\mu}{mu} and variance
    \eqn{\sigma^2}{sigma^2} of \eqn{S} are computed
    from the values of the covariate at all locations in the window.
    Then we compute \eqn{Z_1 = (S-\mu)/\sigma}{Z[1]=(S-mu)/sigma}.
    \item If \code{which="Z2"},
    the test statistic \eqn{Z_2}{Z[2]} is computed as follows.
    The values of the \code{covariate} at all locations in the
    observation window, weighted by the point process intensity,
    are compiled into a cumulative distribution function \eqn{F}.
    The probability integral transformation is then applied:
    the values of the \code{covariate} at the original data points
    are transformed by the predicted cumulative distribution function
    \eqn{F} into numbers between 0 and 1. If the model is correct,
    these numbers are i.i.d. uniform random numbers.
    The standardised sample mean of these numbers is the
    statistic \eqn{Z_2}{Z[2]}. 
  }
  In both cases the null distribution of the test statistic
  is the standard normal distribution, approximately.

  The return value is an object of class \code{"htest"} containing the
  results of the hypothesis test. The print method for this class
  gives an informative summary of the test outcome.
}
\value{
  An object of class \code{"htest"} (hypothesis test)
  and also of class \code{"bermantest"},
  containing the results of the test. The return value can be
  plotted (by \code{\link{plot.bermantest}}) or printed
  to give an informative summary of the test.
}
\section{Warning}{
  The meaning of a one-sided test must be carefully scrutinised: see
  the printed output.
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\seealso{
  \code{\link{kstest}},
  \code{\link{quadrat.test}},
  \code{\link{ppm}}
}
\references{
  Berman, M. (1986)
  Testing for spatial association between a point process
  and another stochastic process.
  \emph{Applied Statistics} \bold{35}, 54--62.
}
\examples{
   # Berman's data
   data(copper)
   X <- copper$SouthPoints
   L <- copper$SouthLines
   D <- distmap(L, eps=1)
   # test of CSR
   bermantest(X, D)
   bermantest(X, D, "Z2")
}
\keyword{htest}
\keyword{spatial}

