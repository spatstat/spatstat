\name{quadrat.test.splitppp}
\alias{quadrat.test.splitppp}
\title{Dispersion Test of CSR for Split Point Pattern Based on
       Quadrat Counts}
\description{
  Performs a test of Complete Spatial Randomness
  for each of the component patterns in a split point pattern,
  based on quadrat counts.
  By default performs chi-squared tests; can also perform
  Monte Carlo based tests.
}
\usage{
\method{quadrat.test}{splitppp}(X, ..., df=NULL, df.est=NULL, Xname=NULL)
}
\arguments{
  \item{X}{
    A split point pattern (object of class \code{"splitppp"}),
    each component of which will be subjected to the goodness-of-fit test.
  }
  \item{\dots}{Arguments passed to \code{\link{quadrat.test.ppp}}.}
  \item{df,df.est,Xname}{Arguments passed to \code{\link{pool.quadrattest}}.}
}
\details{
  The function \code{quadrat.test} is generic, with methods for
  point patterns (class \code{"ppp"}), split point patterns
  (class \code{"splitppp"}) and point process models
  (class \code{"ppm"}).

  If \code{X} is a split point pattern, then for each of the
  component point patterns (taken separately) we test 
  the null hypotheses of Complete Spatial Randomness,
  then combine the result into a single test.

  The method \code{quadrat.test.ppp} is applied to each
  component point pattern. Then the results are pooled using
  \code{\link{pool.quadrattest}} to obtain a single test.
}
\seealso{
  \code{\link{quadrat.test}},
  \code{\link{quadratcount}},
  \code{\link{quadrats}},
  \code{\link{quadratresample}},
  \code{\link{chisq.test}},
  \code{\link{kstest}}.

  To test a Poisson point process model against a specific Poisson alternative,
  use \code{\link{anova.ppm}}.
}
\value{
  An object of class \code{"quadrattest"} which can be printed and
  plotted.
}
\examples{
 data(humberside)
 qH <- quadrat.test(split(humberside), 2, 3)
 plot(qH)
 qH
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{htest}
