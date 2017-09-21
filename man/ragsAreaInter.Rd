\name{ragsAreaInter}
\alias{ragsAreaInter}
\title{
  Alternating Gibbs Sampler for Area-Interaction Process
}
\description{
 Generate a realisation of the area-interaction process
 using the alternating Gibbs sampler.
 Applies only when the interaction parameter \eqn{eta} is greater than 1.
}
\usage{
     ragsAreaInter(beta, eta, r, \dots,
                   win = NULL, bmax = NULL, periodic = FALSE, ncycles = 100)
}
\arguments{
  \item{beta}{
    First order trend. A number, a pixel image (object of class
    \code{"im"}), or a \code{function(x,y)}.
  }
  \item{eta}{
    Interaction parameter (canonical form) as described in
    the help for \code{\link{AreaInter}}.
    A number greater than 1.
  }
  \item{r}{
    Disc radius in the model. A number greater than 1.
  }
  \item{\dots}{
    Additional arguments for \code{beta} if it is a function.
  }
  \item{win}{
    Simulation window. An object of class \code{"owin"}.
    (Ignored if \code{beta} is a pixel image.)
  }
  \item{bmax}{
    Optional. The maximum possible value of \code{beta},
    or a number larger than this. 
  }
  \item{periodic}{
    Logical value indicating whether to treat opposite sides of the
    simulation window as being the same, so that points close to one
    side may interact with points close to the opposite side.
    Feasible only when the window is a rectangle.
  }
  \item{ncycles}{
    Number of cycles of the alternating Gibbs sampler to be performed.
  }
}
\details{
  This function generates a simulated realisation of the
  area-interaction process (see \code{\link{AreaInter}})
  using the alternating Gibbs sampler (see \code{\link{rags}}).

  It exploits a mathematical relationship between the
  (unmarked) area-interaction process and the two-type
  hard core process (Baddeley and Van Lieshout, 1995;
  Widom and Rowlinson, 1970). This relationship only holds
  when the interaction parameter \code{eta} is greater than 1
  so that the area-interaction process is clustered.

  The parameters \code{beta,eta} are the canonical parameters described
  in the help for \code{\link{AreaInter}}.
  The first order trend \code{beta} may be a constant, a function,
  or a pixel image.

  The simulation window is determined by \code{beta} if it is a pixel
  image, and otherwise by the argument \code{win} (the default is the
  unit square).
}
\value{
  A point pattern (object of class \code{"ppp"}).
}
\references{
  Baddeley, A.J. and Van Lieshout, M.N.M. (1995).
  Area-interaction point processes.
  \emph{Annals of the Institute of Statistical Mathematics}
  \bold{47} (1995) 601--619.

  Widom, B. and Rowlinson, J.S. (1970).
  New model for the study of liquid-vapor phase transitions.
  \emph{The Journal of Chemical Physics}
  \bold{52} (1970) 1670--1684.
}
\author{
  \adrian.
}
\seealso{
  \code{\link{rags}},
  \code{\link{ragsMultiHard}}

  \code{\link{AreaInter}}
}
\examples{
   plot(ragsAreaInter(100, 2, 0.07, ncycles=15))
}
\keyword{spatial}
\keyword{datagen}

