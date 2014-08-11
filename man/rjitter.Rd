\name{rjitter}
\alias{rjitter}
\title{Random Perturbation of a Point Pattern}
\description{
  Applies independent random displacements to each point in a point pattern.
}
\usage{
rjitter(X, radius, retry=TRUE, giveup = 10000)
}
\arguments{
  \item{X}{
    A point pattern (object of class \code{"ppp"}).
  }
  \item{radius}{
    Scale of perturbations. A positive numerical value.
    The displacement vectors will be uniformly
    distributed in a circle of this radius.
  }
  \item{retry}{
    What to do when a perturbed point lies outside the window
    of the original point pattern. If \code{retry=FALSE},
    the point will be lost; if \code{retry=TRUE},
    the algorithm will try again.
  }
  \item{giveup}{
    Maximum number of unsuccessful attempts.
  }
}
\details{
  Each of the points in the point pattern \code{X} is subjected to
  an independent random displacement. The displacement vectors are
  uniformly distributed in a circle of radius \code{radius}.

  If a displaced point lies outside the window, then if
  \code{retry=FALSE} the point will be lost. 

  However if \code{retry=TRUE}, the algorithm will try again: each time a
  perturbed point lies outside the window, the algorithm will reject it and
  generate another proposed perturbation of the original point,
  until one lies inside the window, or until \code{giveup} unsuccessful
  attempts have been made. In the latter case, any unresolved points
  will be included without any perturbation. The return value will
  always be a point pattern with the same number of points as \code{X}.
}
\value{
  A point pattern (object of class \code{"ppp"}) in the same window
  as \code{X}. 
}
\examples{
   X <- rsyst(owin(), 10, 10)
   Y <- rjitter(X, 0.02)
   plot(Y)
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{datagen}
