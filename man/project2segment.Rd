\name{project2segment}
\alias{project2segment}
\title{Move Point To Nearest Line}
\description{
  Given a point pattern and a line segment pattern,
  this function moves each point to the closest location
  on a line segment.
}
\usage{
project2segment(X, Y)
}
\arguments{
  \item{X}{A point pattern (object of class \code{"ppp"}).}
  \item{Y}{A line segment pattern (object of class \code{"psp"}).}
}
\details{
  For each point \code{x} in the point pattern \code{X}, this function
  finds the closest line segment \code{y} in the line segment pattern
  \code{Y}. It then `projects' the point \code{x} onto the line segment
  \code{y} by finding the position \code{z} along \code{y}
  which is closest to \code{x}. This position \code{z} is returned,
  along with supplementary information.
}
\value{
  A list with the following components. Each component has length equal to
  the number of points in \code{X}, and its entries correspond to the
  points of \code{X}.
  
  \item{Xproj }{
    Point pattern (object of class \code{"ppp"}
    containing the projected points.
  }
  \item{mapXY }{
    Integer vector identifying the nearest segment to each point.
  }
  \item{d}{
    Numeric vector of distances from each point of \code{X} to
    the corresponding projected point.
  }
  \item{tp}{
    Numeric vector giving the scaled parametric coordinate 
    \eqn{0 \le t_p \le 1}{0 <= tp <= 1} of the position of the projected point
    along the segment.
  }

  For example suppose \code{mapXY[2] = 5} and \code{tp[2] = 0.33}.
  Then \code{Y[5]} is the line segment lying closest to \code{X[2]}.
  The projection of the point \code{X[2]} onto the segment \code{Y[5]}
  is the point \code{Xproj[2]}, which lies one-third of the way
  between the first and second endpoints of the line segment \code{Y[5]}.

  }
\author{
  Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  and Rolf Turner \email{r.turner@auckland.ac.nz}
  
}
\seealso{
  \code{\link{nearestsegment}} for a faster way to determine
  which segment is closest to each point.
}
\examples{
  X <- rstrat(square(1), 5)
  Y <- as.psp(matrix(runif(20), 5, 4), window=owin())
  plot(Y, lwd=3, col="green")
  plot(X, add=TRUE, col="red", pch=16)
  v <- project2segment(X,Y)
  Xproj <- v$Xproj
  plot(Xproj, add=TRUE, pch=16)
  arrows(X$x, X$y, Xproj$x, Xproj$y, angle=10, length=0.15, col="red")
}
\keyword{spatial}
\keyword{math}
