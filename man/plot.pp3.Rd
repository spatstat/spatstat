\name{plot.pp3}
\Rdversion{1.1}
\alias{plot.pp3}
\title{
  Plot a Three-Dimensional Point Pattern
}
\description{
  Plots a three-dimensional point pattern.
}
\usage{
 \method{plot}{pp3}(x, ..., eye=NULL, org=NULL, theta=25, phi=15,
                    type=c("p", "n", "h"))
}
\arguments{
  \item{x}{
    Three-dimensional point pattern (object of class \code{"pp3"}).
  }
  \item{\dots}{
    Arguments passed to \code{\link[graphics]{points}}
    controlling the appearance of the points.
  }
  \item{eye}{
    Optional. Eye position. A numeric vector of length 3 giving the
    location from which the scene is viewed.
  }
  \item{org}{
    Optional. Origin (centre) of the view. A numeric vector of length
    3 which will be at the centre of the view.
  }
  \item{theta,phi}{
    Optional angular coordinates (in degrees) specifying the direction
    from which the scene is viewed: \code{theta} is the azimuth
    and \code{phi} is the colatitude. Ignored if \code{eye} is given.
  }
  \item{type}{
    Type of plot: \code{type="p"} for points,
    \code{type="h"} for points on vertical lines,
    \code{type="n"} for box only.
  }
}
\details{
  This is the plot method for objects of class \code{"pp3"}.
  It generates a two-dimensional plot of the point pattern \code{x}
  and its containing box as if they had been viewed from the
  location specified by \code{eye} (or from the direction
  specified by \code{theta} and \code{phi}). 

  Note that default values of arguments to \code{plot.pp3}
  can be set by \code{\link{spatstat.options}("par.pp3")}.
}
\value{Null.}
\author{
Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}


Rolf Turner \email{r.turner@auckland.ac.nz}

and Ege Rubak \email{rubak@math.aau.dk}

}
\seealso{
  \code{\link{pp3}},
  \code{\link{spatstat.options}}.
}
\examples{
  X <- osteo$pts[[1]]
  plot(X, main="Osteocyte lacunae, animal 1, brick 1",
       cex=1.5, pch=16)
}
\keyword{spatial}
\keyword{hplot}
