\name{rMaternI}
\alias{rMaternI}
\title{Simulate Matern Model I}
\description{
  Generate a random point pattern, a simulated realisation of the
  \ifelse{latex}{\out{Mat\'ern}}{Matern} Model I inhibition process model.
}
\usage{
 rMaternI(kappa, r, win = owin(c(0,1),c(0,1)), stationary=TRUE, \dots,
          nsim=1, drop=TRUE)
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
    Alternatively a higher-dimensional box of class
    \code{"box3"} or \code{"boxx"}.
  }
  \item{stationary}{
    Logical. Whether to start with a stationary process of proposal points
    (\code{stationary=TRUE}) or to generate the
    proposal points only inside the window (\code{stationary=FALSE}).
  }
  \item{\dots}{Ignored.}
  \item{nsim}{Number of simulated realisations to be generated.}
  \item{drop}{
    Logical. If \code{nsim=1} and \code{drop=TRUE} (the default), the
    result will be a point pattern, rather than a list 
    containing a point pattern.
  }
}
\value{
  A point pattern
  if \code{nsim=1}, or a list of point patterns if \code{nsim > 1}.
  Each point pattern is normally an object of class \code{"ppp"},
  but may be of class \code{"pp3"} or \code{"ppx"} depending on the window.
}
\details{
  This algorithm generates one or more realisations
  of \ifelse{latex}{\out{Mat\'ern}}{Matern}'s Model I
  inhibition process inside the window \code{win}.

  The process is constructed by first
  generating a uniform Poisson point process of ``proposal'' points 
  with intensity \code{kappa}. If \code{stationary = TRUE} (the
  default), the proposal points are generated in a window larger than
  \code{win} that effectively means the proposals are stationary.
  If \code{stationary=FALSE} then the proposal points are
  only generated inside the window \code{win}.

  A proposal point is then deleted if it lies within \code{r} units' distance
  of another proposal point. Otherwise it is retained.
  
  The retained points constitute \ifelse{latex}{\out{Mat\'ern}}{Matern}'s Model I.
}
\seealso{
\code{\link{rpoispp}},
\code{\link{rMatClust}}
}
\examples{
 X <- rMaternI(20, 0.05)
 Y <- rMaternI(20, 0.05, stationary=FALSE)
}
\author{
  \adrian
  
  ,
  Ute Hahn,
  \rolf
  
  and \ege
  
}
\keyword{spatial}
\keyword{datagen}

