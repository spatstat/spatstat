\name{rpoint}
\alias{rpoint}
\title{Generate N Random Points}
\description{
  Generate a random point pattern
  containing \eqn{n} independent, identically distributed random points
  with any specified distribution.
}
\usage{
 rpoint(n, f, fmax=NULL, win=unit.square(),
        \dots, giveup=1000, verbose=FALSE,
        nsim=1, drop=TRUE)
}
\arguments{
  \item{n}{
    Number of points to generate.
  }
  \item{f}{
    The probability density of the points,
    possibly un-normalised.
    Either a constant,
    a function \code{f(x,y,...)}, or a pixel image object.
  }
  \item{fmax}{
    An upper bound on the values of \code{f}.
    If missing, this number will be estimated.
  }
  \item{win}{
    Window in which to simulate the pattern.
    Ignored if \code{f} is a pixel image.
  }
  \item{\dots}{
    Arguments passed to the function \code{f}.
  }
  \item{giveup}{
    Number of attempts in the rejection method after which
    the algorithm should stop trying to generate new points.
  }
  \item{verbose}{
    Flag indicating whether to report details of performance
    of the simulation algorithm.
  }
  \item{nsim}{Number of simulated realisations to be generated.}
  \item{drop}{
    Logical. If \code{nsim=1} and \code{drop=TRUE} (the default), the
    result will be a point pattern, rather than a list 
    containing a point pattern.
  }
}
\value{
  A point pattern (an object of class \code{"ppp"})
  if \code{nsim=1}, or a list of point patterns if \code{nsim > 1}.
}
\details{
  This function generates \code{n} independent, identically distributed
  random points with common probability density proportional to
  \code{f}.

  The argument \code{f} may be
  \describe{
    \item{a numerical constant:}{
      uniformly distributed random points will
      be generated.
    }
    \item{a function:}{random points will be generated
      in the window \code{win} with probability density proportional
      to \code{f(x,y,...)} where \code{x} and \code{y} are the cartesian
      coordinates. The function \code{f} must accept 
      two \emph{vectors} of coordinates \code{x,y} and return the corresponding
      vector of function values. Additional arguments \code{...} of any kind
      may be passed to the function.
    }
    \item{a pixel image:}{if \code{f} is a pixel image object
      of class \code{"im"} (see \code{\link{im.object}}) then
      random points will be generated
      in the window of this pixel image, with probability density
      proportional to the pixel values of \code{f}.
    }
  }
  The algorithm is as follows:
  \itemize{
    \item
    If \code{f} is a constant, we invoke \code{\link{runifpoint}}.
    \item
    If \code{f} is a function, then we use the rejection method.
    Proposal points are generated from the uniform distribution.
    A proposal point \eqn{(x,y)} is accepted with probability
    \code{f(x,y,...)/fmax} and otherwise rejected.
    The algorithm continues until \code{n} points have been
    accepted. It gives up after \code{giveup * n} proposals
    if there are still fewer than \code{n} points.
    \item 
    If \code{f} is a pixel image, then a random sequence of 
    pixels is selected (using \code{\link{sample}})
    with probabilities proportional to the
    pixel values of \code{f}.  Then for each pixel in the sequence
    we generate a uniformly distributed random point in that pixel.
  }
  The algorithm for pixel images is more efficient than that for
  functions.
}
\seealso{
\code{\link{ppp.object}},
\code{\link{owin.object}},
\code{\link{runifpoint}}
}
\examples{
 # 100 uniform random points in the unit square
 X <- rpoint(100)

 # 100 random points with probability density proportional to x^2 + y^2
 X <- rpoint(100, function(x,y) { x^2 + y^2}, 1)

 # `fmax' may be omitted
 X <- rpoint(100, function(x,y) { x^2 + y^2})

 # irregular window
 data(letterR)
 X <- rpoint(100, function(x,y) { x^2 + y^2}, win=letterR)

 # make a pixel image 
 Z <- setcov(letterR)
 # 100 points with density proportional to pixel values
 X <- rpoint(100, Z)
}
\author{Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  and Rolf Turner \email{r.turner@auckland.ac.nz}
  
}
\keyword{spatial}
\keyword{datagen}
