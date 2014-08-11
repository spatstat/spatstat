\name{rMaternII}
\alias{rMaternII}
\title{Simulate Matern Model II}
\description{
  Generate a random point pattern, a simulated realisation of the
  \ifelse{latex}{\out{Mat\'ern}}{Matern} Model II inhibition process.
}
\usage{
 rMaternII(kappa, r, win = owin(c(0,1),c(0,1)), stationary=TRUE)
}
\arguments{
  \item{kappa}{
    Intensity of the Poisson process of proposal points.
    A single positive number.
  }
  \item{r}{
    Inhibition distance.
  }
  \item{win}{
    Window in which to simulate the pattern.
    An object of class \code{"owin"}
    or something acceptable to \code{\link{as.owin}}.
  }
  \item{stationary}{
    Logical. Whether to start with a stationary process of proposal points
    (\code{stationary=TRUE}) or to generate the
    proposal points only inside the window (\code{stationary=FALSE}).
  }
}
\value{
  The simulated point pattern (an object of class \code{"ppp"}).
}
\details{
  This algorithm generates a realisation of \ifelse{latex}{\out{Mat\'ern}}{Matern}'s Model II
  inhibition process 
  inside the window \code{win}.

  The process is constructed by first
  generating a uniform Poisson point process of ``proposal'' points 
  with intensity \code{kappa}. If \code{stationary = TRUE} (the
  default), the proposal points are generated in a window larger than
  \code{win} that effectively means the proposals are stationary.
  If \code{stationary=FALSE} then the proposal points are
  only generated inside the window \code{win}.
  
  Then each proposal point is marked by an ``arrival time'', a number
  uniformly distributed in \eqn{[0,1]} independently of other variables.
  
  A proposal point is deleted if it lies within \code{r} units' distance
  of another proposal point \emph{that has an earlier arrival time}.
  Otherwise it is retained.
  The retained points constitute \ifelse{latex}{\out{Mat\'ern}}{Matern}'s Model II.

  The difference between \ifelse{latex}{\out{Mat\'ern}}{Matern}'s Model I and II is the italicised
  statement above. Model II has a higher intensity
  for the same parameter values.
}
\seealso{
\code{\link{rpoispp}},
\code{\link{rMatClust}},
\code{\link{rMaternI}}
}
\examples{
 X <- rMaternII(20, 0.05)
 Y <- rMaternII(20, 0.05, stationary=FALSE)
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/},
  Ute Hahn, 
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{datagen}

