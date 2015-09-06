\name{rmpoispp}
\alias{rmpoispp}
\title{Generate Multitype Poisson Point Pattern}
\description{
  Generate a random point pattern, a realisation of the
  (homogeneous or inhomogeneous) multitype Poisson process.
}
\usage{
 rmpoispp(lambda, lmax=NULL, win, types, \dots,
          nsim=1, drop=TRUE, warnwin=!missing(win))
}
\arguments{
  \item{lambda}{
    Intensity of the multitype Poisson process.
    Either a single positive number, a vector, a \code{function(x,y,m, \dots)},
    a pixel image, a list of functions \code{function(x,y, \dots)},
    or a list of pixel images.
  }
  \item{lmax}{
    An upper bound for the value of \code{lambda}.
    May be omitted
  }
  \item{win}{
    Window in which to simulate the pattern.
    An object of class \code{"owin"}
    or something acceptable to \code{\link{as.owin}}.
    Ignored if \code{lambda} is a pixel image or list of images.
  }
  \item{types}{
    All the possible types for the multitype pattern.
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
  \item{warnwin}{
    Logical value specifying whether to issue a warning
    when \code{win} is ignored.
  }
}
\value{
  A point pattern (an object of class \code{"ppp"}) if \code{nsim=1},
  or a list of point patterns if \code{nsim > 1}.
  Each point pattern is multitype (it carries a vector of marks
  which is a factor).
}
\details{
  This function generates a realisation of the marked Poisson
  point process with intensity \code{lambda}.
  
  Note that the intensity function
  \eqn{\lambda(x,y,m)}{lambda(x,y,m)} is the
  average number of points \bold{of type m} per unit area
  near the location \eqn{(x,y)}.
  Thus a marked point process with a constant intensity of 10
  and three possible types will have an average of 30 points per unit
  area, with 10 points of each type on average.

  The intensity function may be specified in any of the following ways.
  \describe{
    \item{single number:}{
      If \code{lambda} is a single number,
      then this algorithm generates a realisation
      of the uniform marked Poisson process inside the window \code{win} with 
      intensity \code{lambda} for each type. The total intensity of
      points of all types is \code{lambda * length(types)}.
      The argument \code{types} must be given
      and determines the possible types in the multitype pattern.
    }
    \item{vector:}{
      If \code{lambda} is a numeric vector,
      then this algorithm generates a realisation
      of the stationary marked Poisson process inside the window
      \code{win} with intensity \code{lambda[i]} for points of type
      \code{types[i]}. The total intensity of points of all types
      is \code{sum(lambda)}.
      The argument \code{types} defaults to
      \code{names(lambda)}, or if that is null, \code{1:length(lambda)}.
    }
    \item{function:}{
      If \code{lambda} is a function, the process has intensity
      \code{lambda(x,y,m,\dots)} at spatial location \code{(x,y)}
      for points of type \code{m}.
      The function \code{lambda} must work correctly with vectors \code{x},
      \code{y} and \code{m}, returning a vector of function values.
      (Note that \code{m} will be a factor
      with levels equal to \code{types}.)
      The value \code{lmax}, if present, must be an upper bound on the 
      values of \code{lambda(x,y,m,\dots)} for all locations \code{(x, y)}
      inside the window \code{win} and all types \code{m}.
      The argument \code{types} must be given.
    }
    \item{list of functions:}{
      If \code{lambda} is a list of functions,
      the process has intensity \code{lambda[[i]](x,y,\dots)} at spatial
      location \code{(x,y)} for points of type \code{types[i]}.
      The function \code{lambda[[i]]} must work correctly with vectors
      \code{x} and \code{y}, returning a vector of function values.
      The value \code{lmax}, if given, must be an upper bound on the 
      values of \code{lambda(x,y,\dots)} for all locations \code{(x, y)}
      inside the window \code{win}. 
      The argument \code{types} defaults to
      \code{names(lambda)}, or if that is null, \code{1:length(lambda)}.
    }
    \item{pixel image:}{
      If \code{lambda} is a pixel image object of class \code{"im"}
      (see \code{\link{im.object}}), the intensity at a location
      \code{(x,y)} for points of any type is equal to the pixel value
      of \code{lambda} for the pixel nearest to \code{(x,y)}.
      The argument \code{win} is ignored;
      the window of the pixel image is used instead.
      The argument \code{types} must be given.
    }
    \item{list of pixel images:}{
      If \code{lambda} is a list of pixel images,
      then the image \code{lambda[[i]]} determines the intensity
      of points of type \code{types[i]}. 
      The argument \code{win} is ignored;
      the window of the pixel image is used instead.
      The argument \code{types} defaults to
      \code{names(lambda)}, or if that is null, \code{1:length(lambda)}.
    }
  }
  If \code{lmax} is missing, an approximate upper bound will be calculated.
  
  To generate an inhomogeneous Poisson process
  the algorithm uses ``thinning'': it first generates a uniform
  Poisson process of intensity \code{lmax} for points of each type \code{m},
  then randomly deletes or retains each point independently,
  with retention probability
  \eqn{p(x,y,m) = \lambda(x,y,m)/\mbox{lmax}}{p(x,y,m) = lambda(x,y)/lmax}.
}
\seealso{
  \code{\link{rpoispp}} for unmarked Poisson point process;
  \code{\link{rmpoint}} for a fixed number of random marked points;
  \code{\link{ppp.object}},
  \code{\link{owin.object}}.
}
\examples{
 # uniform bivariate Poisson process with total intensity 100 in unit square
 pp <- rmpoispp(50, types=c("a","b"))
 
 # stationary bivariate Poisson process with intensity A = 30, B = 70
 pp <- rmpoispp(c(30,70), types=c("A","B"))
 pp <- rmpoispp(c(30,70))

 # works in any window
 data(letterR)
 pp <- rmpoispp(c(30,70), win=letterR, types=c("A","B"))

 # inhomogeneous lambda(x,y,m)
 # note argument 'm' is a factor 
 lam <- function(x,y,m) { 50 * (x^2 + y^3) * ifelse(m=="A", 2, 1)}
 pp <- rmpoispp(lam, win=letterR, types=c("A","B"))
 # extra arguments
 lam <- function(x,y,m,scal) { scal * (x^2 + y^3) * ifelse(m=="A", 2, 1)}
 pp <- rmpoispp(lam, win=letterR, types=c("A","B"), scal=50)

 # list of functions lambda[[i]](x,y)
 lams <- list(function(x,y){50 * x^2}, function(x,y){20 * abs(y)})
 pp <- rmpoispp(lams, win=letterR, types=c("A","B"))
 pp <- rmpoispp(lams, win=letterR)
 # functions with extra arguments
 lams <- list(function(x,y,scal){5 * scal * x^2},
              function(x,y, scal){2 * scal * abs(y)})
 pp <- rmpoispp(lams, win=letterR, types=c("A","B"), scal=10)
 pp <- rmpoispp(lams, win=letterR, scal=10)

 # florid example
 lams <- list(function(x,y){
   			   100*exp((6*x + 5*y - 18*x^2 + 12*x*y - 9*y^2)/6)
                         }
                         # log quadratic trend
              ,
              function(x,y){
                         	   100*exp(-0.6*x+0.5*y)
                         }
                        # log linear trend
              )
  X <- rmpoispp(lams, win=unit.square(), types=c("on", "off"))   

 # pixel image
 Z <- as.im(function(x,y){30 * (x^2 + y^3)}, letterR)
 pp <- rmpoispp(Z, types=c("A","B"))

 # list of pixel images
 ZZ <- list(
          as.im(function(x,y){20 * (x^2 + y^3)}, letterR),
          as.im(function(x,y){40 * (x^3 + y^2)}, letterR))
 pp <- rmpoispp(ZZ, types=c("A","B"))
 pp <- rmpoispp(ZZ)

 # randomising an existing point pattern
 rmpoispp(intensity(amacrine), win=Window(amacrine))
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{datagen}

