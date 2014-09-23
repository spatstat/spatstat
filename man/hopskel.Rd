\name{hopskel}
\alias{hopskel}
\alias{hopskel.test}
\title{Hopkins-Skellam Test}
\description{
  Perform the Hopkins-Skellam test of Complete Spatial Randomness,
  or simply calculate the test statistic.
}
\usage{
hopskel(X)

hopskel.test(X, \dots,
             alternative=c("two.sided", "less", "greater",
                           "clustered", "regular"),
             method=c("asymptotic", "MonteCarlo"),
             nsim=999)
}
\arguments{
  \item{X}{
    Point pattern (object of class \code{"ppp"}).
  }
  \item{alternative}{
    String indicating the type of alternative for the
    hypothesis test. Partially matched.
  }
  \item{method}{
    Method of performing the test. Partially matched.
  }
  \item{nsim}{
    Number of Monte Carlo simulations to perform, if a Monte Carlo
    p-value is required.
  }
  \item{\dots}{Ignored.}
}
\details{
  Hopkins and Skellam (1954) proposed a test of Complete Spatial
  Randomness based on comparing nearest-neighbour distances with
  point-event distances.

  If the point pattern \code{X} contains \code{n}
  points, we first compute the nearest-neighbour distances
  \eqn{P_1, \ldots, P_n}{P[1], ..., P[n]} 
  so that \eqn{P_i}{P[i]} is the distance from the \eqn{i}th data
  point to the nearest other data point. Then we 
  generate another completely random pattern \code{U} with
  the same number \code{n} of points, and compute for each point of \code{U}
  the distance to the nearest point of \code{X}, giving
  distances \eqn{I_1, \ldots, I_n}{I[1], ..., I[n]}.
  The test statistic is 
  \deqn{
    A = \frac{\sum_i P_i^2}{\sum_i I_i^2}
  }{
    A = (sum[i] P[i]^2) / (sum[i] I[i]^2)
  }
  The null distribution of \eqn{A} is roughly
  an \eqn{F} distribution with shape parameters \eqn{(2n,2n)}.
  (This is equivalent to using the test statistic \eqn{H=A/(1+A)}
  and referring \eqn{H} to the Beta distribution with parameters
  \eqn{(n,n)}).

  The function \code{hopskel} calculates the Hopkins-Skellam test statistic
  \eqn{A}, and returns its numeric value. This can be used as a simple
  summary of spatial pattern: the value \eqn{H=1} is consistent
  with Complete Spatial Randomness, while values \eqn{H < 1} are
  consistent with spatial clustering, and values \eqn{H > 1} are consistent
  with spatial regularity.

  The function \code{hopskel.test} performs the test.
  If \code{method="asymptotic"} (the default), the test statistic \eqn{H}
  is referred to the \eqn{F} distribution. If \code{method="MonteCarlo"},
  a Monte Carlo test is performed using \code{nsim} simulated point
  patterns.
}
\value{
  The value of \code{hopskel} is a single number.

  The value of \code{hopskel.test} is an object of class \code{"htest"}
  representing the outcome of the test. It can be printed. 
}
\references{
  Hopkins, B. and Skellam, J.G. (1954) 
  A new method of determining the type of distribution
  of plant individuals. \emph{Annals of Botany} \bold{18}, 
  213--227.
}
\seealso{
  \code{\link{clarkevans}},
  \code{\link{clarkevans.test}},
  \code{\link{nndist}},
  \code{\link{nncross}}
}
\examples{
  hopskel(redwood)
  hopskel(redwood)
  hopskel.test(redwood, alternative="clustered")
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  Rolf Turner
  \email{r.turner@auckland.ac.nz}
  and Ege Rubak
  \email{rubak@math.aau.dk}
}
\keyword{spatial}
\keyword{nonparametric}
\keyword{htest}
