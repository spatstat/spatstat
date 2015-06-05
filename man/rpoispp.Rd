\name{rpoispp}
\alias{rpoispp}
\title{Generate Poisson Point Pattern}
\description{
  Generate a random point pattern using the
  (homogeneous or inhomogeneous) Poisson process.
  Includes CSR (complete spatial randomness).
}
\usage{
 rpoispp(lambda, lmax=NULL, win=owin(), \dots, nsim=1, drop=TRUE, ex=NULL)
}
\arguments{
  \item{lambda}{
    Intensity of the Poisson process.
    Either a single positive number, a \code{function(x,y, \dots)},
    or a pixel image.
  }
  \item{lmax}{
    An upper bound for the value of \code{lambda(x,y)},
    if \code{lambda} is a function.
  }
  \item{win}{
    Window in which to simulate the pattern.
    An object of class \code{"owin"}
    or something acceptable to \code{\link{as.owin}}.
    Ignored if \code{lambda} is a pixel image.
  }
  \item{\dots}{
    Arguments passed to \code{lambda} if it is a function.
  }
  \item{nsim}{Number of simulated realisations to be generated.}
  \item{drop}{
    Logical. If \code{nsim=1} and \code{drop=TRUE} (the default), the
    result will be a point pattern, rather than a list 
    containing a point pattern.
  }
  \item{ex}{
    Optional. A point pattern to use as the example.
    If \code{ex} is given and \code{lambda,lmax,win} are missing,
    then \code{lambda} and \code{win} will be calculated from
    the point pattern \code{ex}.
  }
}
\value{
  A point pattern (an object of class \code{"ppp"})
  if \code{nsim=1}, or a list of point patterns if \code{nsim > 1}.
}
\details{
  If \code{lambda} is a single number,
  then this algorithm generates a realisation
  of the uniform Poisson process (also known as 
  Complete Spatial Randomness, CSR) inside the window \code{win} with 
  intensity \code{lambda} (points per unit area).
 
  If \code{lambda} is a function, then this algorithm generates a realisation
  of the inhomogeneous Poisson process with intensity function
  \code{lambda(x,y,\dots)} at spatial location \code{(x,y)}
  inside the window \code{win}.
  The function \code{lambda} must work correctly with vectors \code{x}
  and \code{y}.
  
  If \code{lmax} is given,
  it must be an upper bound on the values of \code{lambda(x,y,\dots)}
  for all locations \code{(x, y)}
  inside the window \code{win}.

  If \code{lmax} is missing or \code{NULL},
  an approximate upper bound is computed by finding the maximum value
  of \code{lambda(x,y,\dots)}
  on a grid of locations \code{(x,y)} inside the window \code{win},
  and adding a safety margin equal to 5 percent of the range of
  \code{lambda} values. This can be computationally intensive,
  so it is advisable to specify \code{lmax} if possible.

  If \code{lambda} is a pixel image object of class \code{"im"}
  (see \code{\link{im.object}}), this algorithm generates a realisation
  of the inhomogeneous Poisson process with intensity equal to the
  pixel values of the image. (The value of the intensity function at an
  arbitrary location is the pixel value of the nearest pixel.)
  The argument \code{win} is ignored;
  the window of the pixel image is used instead. It will be converted
  to a rectangle if possible, using \code{\link{rescue.rectangle}}.
  
  To generate an inhomogeneous Poisson process
  the algorithm uses ``thinning'': it first generates a uniform
  Poisson process of intensity \code{lmax},
  then randomly deletes or retains each point, independently of other points,
  with retention probability
  \eqn{p(x,y) = \lambda(x,y)/\mbox{lmax}}{p(x,y) = lambda(x,y)/lmax}.

  For \emph{marked} point patterns, use \code{\link{rmpoispp}}.
}
\section{Warning}{
  Note that \code{lambda} is the \bold{intensity}, that is,
  the expected number of points \bold{per unit area}.
  The total number of points in the simulated
  pattern will be random with expected value \code{mu = lambda * a}
  where \code{a} is the area of the window \code{win}. 
}
\seealso{
  \code{\link{rmpoispp}} for Poisson \emph{marked} point patterns,
  \code{\link{runifpoint}} for a fixed number of independent
  uniform random points;
  \code{\link{rpoint}}, \code{\link{rmpoint}} for a fixed number of
  independent random points with any distribution;
  \code{\link{rMaternI}},
  \code{\link{rMaternII}},
  \code{\link{rSSI}},
  \code{\link{rStrauss}},
  \code{\link{rstrat}}
  for random point processes with spatial inhibition
  or regularity; 
  \code{\link{rThomas}},
  \code{\link{rGaussPoisson}},
  \code{\link{rMatClust}},
  \code{\link{rcell}}
  for random point processes exhibiting clustering;
  \code{\link{rmh.default}} for Gibbs processes.
  See also \code{\link{ppp.object}},
  \code{\link{owin.object}}.
}
\examples{
 # uniform Poisson process with intensity 100 in the unit square
 pp <- rpoispp(100)
 
 # uniform Poisson process with intensity 1 in a 10 x 10 square
 pp <- rpoispp(1, win=owin(c(0,10),c(0,10)))
 # plots should look similar !
 
 # inhomogeneous Poisson process in unit square
 # with intensity lambda(x,y) = 100 * exp(-3*x)
 # Intensity is bounded by 100
 pp <- rpoispp(function(x,y) {100 * exp(-3*x)}, 100)

 # How to tune the coefficient of x
 lamb <- function(x,y,a) { 100 * exp( - a * x)}
 pp <- rpoispp(lamb, 100, a=3)

 # pixel image
 Z <- as.im(function(x,y){100 * sqrt(x+y)}, unit.square())
 pp <- rpoispp(Z)

 # randomising an existing point pattern
 rpoispp(intensity(cells), win=Window(cells))
 rpoispp(ex=cells)
}
\author{
Adrian Baddeley
\email{Adrian.Baddeley@uwa.edu.au}
\url{http://www.maths.uwa.edu.au/~adrian/}
Rolf Turner
\email{r.turner@auckland.ac.nz}
and Ege Rubak
\email{rubak@math.aau.dk}
}
\keyword{spatial}
\keyword{datagen}
