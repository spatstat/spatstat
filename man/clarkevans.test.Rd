\name{clarkevans.test}
\alias{clarkevans.test}
\title{Clark and Evans Test}
\description{
  Performs the Clark-Evans test of aggregation
  for a spatial point pattern.
}
\usage{
clarkevans.test(X, ...,
               correction="none",
               clipregion=NULL,
               alternative=c("two.sided", "less", "greater"),
               nsim=1000)
}
\arguments{
  \item{X}{
    A spatial point pattern (object of class \code{"ppp"}).
  }
  \item{\dots}{Ignored.}
  \item{correction}{
    Character string.
    The type of edge correction to be applied.
    See \code{\link{clarkevans}}
  }
  \item{clipregion}{
    Clipping region for the guard area correction.
    A window (object of class \code{"owin"}).
    See \code{\link{clarkevans}}
  }
  \item{alternative}{
    String indicating the type of alternative for the
    hypothesis test.
  }
  \item{nsim}{
    Number of Monte Carlo simulations to perform, if a Monte Carlo
    p-value is required.
  }
}
\details{
  This command uses the Clark and Evans (1954) aggregation index \eqn{R}
  as the basis for a crude test of clustering or ordering of a point pattern.
  
  The Clark-Evans index is computed by the function
  \code{\link{clarkevans}}. See the help for \code{\link{clarkevans}}
  for information about the Clark-Evans index \eqn{R} and about
  the arguments \code{correction} and \code{clipregion}.

  This command performs a hypothesis test of clustering or ordering of
  the point pattern \code{X}. The null hypothesis is Complete
  Spatial Randomness, i.e.\ a uniform Poisson process. The alternative
  hypothesis is specified by the argument \code{alternative}:
  \itemize{
    \item \code{alternative="less"} or \code{alternative="clustered"}:
    the alternative hypothesis
    is that \eqn{R < 1} corresponding to a clustered point pattern;
    \item \code{alternative="greater"} or \code{alternative="regular"}:
    the alternative hypothesis
    is that \eqn{R > 1} corresponding to a regular or ordered point pattern;
    \item \code{alternative="two.sided"}:
    the alternative hypothesis is that \eqn{R \neq 1}{R != 1}
    corresponding to a clustered or regular pattern.
  }
  
  The Clark-Evans index \eqn{R} is computed for the data
  as described in \code{\link{clarkevans}}.

  If \code{correction="none"}, 
  the \eqn{p}-value for the test is computed by standardising
  \eqn{R} as proposed by Clark and Evans (1954) and referring the
  statistic to the standard Normal distribution.

  For other edge corrections, the \eqn{p}-value for the test is computed
  by Monte Carlo simulation of \code{nsim} realisations of
  Complete Spatial Randomness.
}
\value{
  An object of class \code{"htest"} representing the result of the test.
}
\references{
  Clark, P.J. and Evans, F.C. (1954)
  Distance to nearest neighbour as a measure of spatial
  relationships in populations \emph{Ecology} \bold{35},
  445--453.
  
  Donnelly, K. (1978) Simulations to determine the variance
  and edge-effect of total nearest neighbour distance.
  In \emph{Simulation methods in archaeology},
  Cambridge University Press, pp 91--95.
}
\author{
  Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
}
\seealso{
  \code{\link{clarkevans}}
}
\examples{
  # Example of a clustered pattern
  data(redwood)
  clarkevans.test(redwood)
  clarkevans.test(redwood, alternative="less")
}
\keyword{spatial}
\keyword{nonparametric}
