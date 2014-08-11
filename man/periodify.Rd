\name{periodify} %DontDeclareMethods
\alias{periodify}
\alias{periodify.owin}
\alias{periodify.ppp}
\alias{periodify.psp}
\title{
  Make Periodic Copies of a Spatial Pattern  
}
\description{
  Given a spatial pattern (point pattern, line segment pattern,
  window, etc) make shifted copies of the pattern
  and optionally combine them to make a periodic pattern.
}
\usage{
periodify(X, ...)
\method{periodify}{ppp}(X, nx = 1, ny = 1, ...,
           combine=TRUE, warn=TRUE, check=TRUE,
           ix=(-nx):nx, iy=(-ny):ny,
           ixy=expand.grid(ix=ix,iy=iy))
\method{periodify}{psp}(X, nx = 1, ny = 1, ..., 
           combine=TRUE, warn=TRUE, check=TRUE,
           ix=(-nx):nx, iy=(-ny):ny, 
           ixy=expand.grid(ix=ix,iy=iy))
\method{periodify}{owin}(X, nx = 1, ny = 1, ...,
           combine=TRUE, warn=TRUE, 
           ix=(-nx):nx, iy=(-ny):ny,
           ixy=expand.grid(ix=ix,iy=iy))
}
\arguments{
  \item{X}{
    An object representing a spatial pattern
    (point pattern, line segment pattern or window).
  }
  \item{nx,ny}{
    Integers.
    Numbers of additional copies of \code{X} in each direction.
    The result will be a grid of \code{2 * nx + 1} by \code{2 * ny + 1}
    copies of the original object.
    (Overruled by \code{ix, iy, ixy}).
  }
  \item{\dots}{
    Ignored.
  }
  \item{combine}{
    Logical flag determining whether the copies should be superimposed
    to make an object like \code{X} (if \code{combine=TRUE}) or
    simply returned as a list of objects (\code{combine=FALSE}).
  }
  \item{warn}{
    Logical flag determining whether to issue warnings.
  }
  \item{check}{
    Logical flag determining whether to check the validity of the
    combined pattern.
  }
  \item{ix, iy}{
    Integer vectors determining the grid positions of the copies
    of \code{X}. (Overruled by \code{ixy}).
  }
  \item{ixy}{
    Matrix or data frame with two columns, giving the 
    grid positions of the copies of \code{X}. 
  }
}
\details{
  Given a spatial pattern (point pattern, line segment pattern, etc)
  this function makes a number of shifted copies of the pattern
  and optionally combines them. The function \code{periodify} is
  generic, with methods for various kinds of spatial objects.

  The default is to make a 3 by 3 array of copies of \code{X} and
  combine them into a single pattern of the same kind as \code{X}.
  This can be used (for example) to compute toroidal or periodic
  edge corrections for various operations on \code{X}.

  If the arguments \code{nx}, \code{ny} are given
  and other arguments are missing,
  the original object will be copied \code{nx} times to the right
  and \code{nx} times to the left, then \code{ny} times upward and
  \code{ny} times downward, making \code{(2 * nx + 1) * (2 * ny + 1)} 
  copies altogether, arranged in a grid, centred on the original object.
  
  If the arguments \code{ix}, \code{iy} or \code{ixy} are specified,
  then these determine the grid positions of the copies of \code{X}
  that will be made. For example \code{(ix,iy) = (1, 2)} means a
  copy of \code{X} shifted by the vector \code{(ix * w, iy * h)} where
  \code{w,h} are the width and height of the bounding rectangle of \code{X}.

  If \code{combine=TRUE} (the default) the copies of \code{X} are
  superimposed to create an object of the same kind as \code{X}.
  If \code{combine=FALSE} the copies of \code{X} are returned as a list.
}
\value{
  If \code{combine=TRUE}, an object of the same class as \code{X}.
  If \code{combine=FALSE}, a list of objects of the same class as \code{X}.
}
\seealso{
  \code{\link{shift}}
}
\examples{
  data(cells)
  plot(periodify(cells))
  a <- lapply(periodify(cells$window, combine=FALSE),
        plot, add=TRUE,lty=2)

}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{manip}
