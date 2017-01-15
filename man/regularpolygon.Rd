\name{regularpolygon}
\alias{regularpolygon}
\alias{hexagon}
\title{
  Create A Regular Polygon
}
\description{
  Create a window object representing a regular (equal-sided) polygon.
}
\usage{
regularpolygon(n, edge = 1, centre = c(0, 0), \dots,
               align = c("bottom", "top", "left", "right", "no"))

hexagon(edge = 1, centre = c(0,0), \dots,
        align = c("bottom", "top", "left", "right", "no"))
}
\arguments{
  \item{n}{
    Number of edges in the polygon.
  }
  \item{edge}{
    Length of each edge in the polygon. A single positive number.
  }
  \item{centre}{
    Coordinates of the centre of the polygon.
    A numeric vector of length 2,
    or a \code{list(x,y)} giving the coordinates of exactly one point, or a
    point pattern (object of class \code{"ppp"}) containing exactly one point.
  }
  \item{align}{
    Character string specifying whether to align one of the edges with
    a vertical or horizontal boundary.
  }
  \item{\dots}{
    Ignored.
  }
}
\details{
  The function \code{regularpolygon} creates a regular (equal-sided)
  polygon with \code{n} sides, centred at \code{centre},
  with sides of equal length \code{edge}.
  The function \code{hexagon} is the special case \code{n=6}.

  The orientation of the polygon is determined by the argument \code{align}.
  If \code{align="no"}, one vertex of the polygon is placed on the
  \eqn{x}-axis.
  Otherwise, an edge of the polygon is aligned with one side of
  the frame, specified by the value of \code{align}.
}
\value{
  A window (object of class \code{"owin"}).
}
\author{
  \spatstatAuthors.
}
\seealso{
  \code{\link{disc}}, \code{\link{ellipse}},
  \code{\link{owin}}.

  \code{\link{hextess}} for hexagonal tessellations.
}
\examples{
  plot(hexagon())
  plot(regularpolygon(7))
  plot(regularpolygon(7, align="left"))
}
\keyword{spatial}
\keyword{datagen}
