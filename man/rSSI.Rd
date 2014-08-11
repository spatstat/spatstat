\name{rSSI}
\alias{rSSI}
\title{Simulate Simple Sequential Inhibition}
\description{
  Generate a random point pattern, a realisation of the
  Simple Sequential Inhibition (SSI) process.
}
\usage{
 rSSI(r, n=Inf, win = square(1), giveup = 1000, x.init=NULL)
}
\arguments{
  \item{r}{
    Inhibition distance.
  }
  \item{n}{
    Maximum number of points allowed. 
    If \code{n} is finite, stop when the \emph{total} number of points
    in the point pattern reaches \code{n}.
    If \code{n} is infinite (the default), stop only when
    it is apparently impossible to add any more points.
    See \bold{Details}.
  }
  \item{win}{
    Window in which to simulate the pattern.
    An object of class \code{"owin"}
    or something acceptable to \code{\link{as.owin}}.
    The default window is the unit square, unless 
    \code{x.init} is specified, when the default window 
    is the window of \code{x.init}.
  }
  \item{giveup}{
    Number of rejected proposals after which the algorithm should terminate.
  }
  \item{x.init}{
    Optional. Initial configuration of points. A point pattern
    (object of class \code{"ppp"}).  The pattern returned by
    \code{rSSI} consists of this pattern together with the points
    added via simple sequential inhibition.  See \bold{Details}.
  }
}
\value{
  The simulated point pattern (an object of class \code{"ppp"}).
}
\details{
  This algorithm generates a realisation of the Simple Sequential
  Inhibition point process inside the window \code{win}.

  Starting with an empty window (or with the point pattern
  \code{x.init} if specified), the algorithm adds points
  one-by-one. Each new point is generated uniformly in the window
  and independently of preceding points. If the new point lies
  closer than \code{r} units from an existing point, then it is
  rejected and another random point is generated.
  The algorithm terminates when either
  \describe{
    \item{(a)}{
      the desired number \code{n} of points is reached, or
    }
    \item{(b)}{
      the current point configuration
      has not changed for \code{giveup} iterations,
      suggesting that it is no longer possible to add new points.
    }
  }
  If \code{n} is infinite (the default) then the algorithm terminates
  only when (b) occurs. The result is sometimes called a
  \emph{Random Sequential Packing}.

  Note that argument \code{n} specifies the maximum permitted
  \bold{total} number of points in the pattern returned by
  \code{rSSI()}. If \code{x.init} is not \code{NULL} then
  the number of points that are \emph{added}
  is at most \code{n - npoints(x.init)} if \code{n} is finite.

  Thus if \code{x.init} is not \code{NULL} then argument \code{n}
  must be at least as large as \code{npoints(x.init)}, otherwise
  an error is given.  If \code{n==npoints(x.init)} then a warning
  is given and the call to \code{rSSI()} has no real effect;
  \code{x.init} is returned.

  There is no requirement that the points of \code{x.init} be at
  a distance at least \code{r} from each other.  All of the \emph{added}
  points will be at a distance at least \code{r} from each other
  and from any point of \code{x.init}.

  The points will be generated inside the window \code{win}
  and the result will be a point pattern in the same window.
  
  The default window is the unit square, \code{win = square(1)},
  unless \code{x.init} is specified, when the default
  is \code{win=as.owin(x.init)}, the window of \code{x.init}.

  If both \code{win} and \code{x.init} are specified, and if the
  two windows are different, then a warning will be issued.
  Any points of \code{x.init} lying outside \code{win} will be removed,
  with a warning.
}
\seealso{
\code{\link{rpoispp}},
\code{\link{rMaternI}},
\code{\link{rMaternII}}.
}
\examples{
 Vinf <- rSSI(0.05)

 V100 <- rSSI(0.05, 100)

 X <- runifpoint(100)
 Y <- rSSI(0.03,142,x.init=X) # Y consists of X together with
                              # 42 added points.
 \dontrun{
  plot(Y)
  plot(X,add=TRUE,chars=20,cols="red")
 }
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{datagen}

