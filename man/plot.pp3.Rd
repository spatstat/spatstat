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
                    type=c("p", "n", "h"),
                    box.back=list(col="pink"),
                    box.front=list(col="blue", lwd=2))
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
  \item{box.front,box.back}{
    How to plot the three-dimensional box that contains the points.
    A list of graphical arguments passed to \code{\link[graphics]{segments}},
    or a logical value indicating whether or not to plot
    the relevant part of the box. See Details.
  }
}
\details{
  This is the plot method for objects of class \code{"pp3"}.
  It generates a two-dimensional plot of the point pattern \code{x}
  and its containing box as if they had been viewed from the
  location specified by \code{eye} (or from the direction
  specified by \code{theta} and \code{phi}). 

  The edges of the box at the \sQuote{back} of the scene (as viewed from the
  eye position) are plotted first. Then the points are added. Finally the
  remaining \sQuote{front} edges are plotted. The arguments
  \code{box.back} and \code{box.front} specify graphical parameters for
  drawing the back and front edges, respectively. Alternatively
  \code{box.back=FALSE} specifies that the back edges shall not be drawn.
  
  Note that default values of arguments to \code{plot.pp3}
  can be set by \code{\link{spatstat.options}("par.pp3")}.
}
\value{Null.}
\author{
  \spatstatAuthors.
}
\seealso{
  \code{\link{pp3}},
  \code{\link{spatstat.options}}.
}
\examples{
  X <- osteo$pts[[1]]
  plot(X, main="Osteocyte lacunae, animal 1, brick 1",
       cex=1.5, pch=16)
  plot(X, type="h", main="", box.back=list(lty=3))
}
\keyword{spatial}
\keyword{hplot}
