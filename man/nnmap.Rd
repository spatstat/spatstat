\name{nnmap}
\alias{nnmap}
\title{
  K-th Nearest Point Map
}
\description{
  Given a point pattern, this function constructs pixel images
  giving the distance from each pixel to its \eqn{k}-th nearest
  neighbour in the point pattern, and the index of the \eqn{k}-th nearest
  neighbour.
}
\usage{
  nnmap(X, k = 1, what = c("dist", "which"),
  \dots, W = as.owin(X),
  is.sorted.X = FALSE, sortby = c("range", "var", "x", "y"))
}
\arguments{
  \item{X}{
    Point pattern (object of class \code{"ppp"}).
  }
  \item{k}{
    Integer, or integer vector. The algorithm will find the 
    \code{k}th nearest neighbour. 
  }
  \item{what}{
    Character string specifying what information should be returned.
    Either the nearest neighbour distance (\code{"dist"}),
    the index of the nearest neighbour (\code{"which"}),
    or both.
  }
  \item{\dots}{
    Arguments passed to \code{\link{as.mask}} to determine the pixel
    resolution of the result.
  }
  \item{W}{
    Window (object of class \code{"owin"}) specifying the spatial domain
    in which the distances will be computed.
    Defaults to the window of \code{X}.
  }
  \item{is.sorted.X}{
    Logical value attesting whether the point pattern \code{X}
    has been sorted. See Details.
  }
  \item{sortby}{
    Determines which coordinate to use to sort the point pattern.
    See Details.
  }
}
\details{
  Given a point pattern \code{X}, this function constructs two pixel images:
  \itemize{
    \item a distance map giving, for each pixel,
    the distance to the nearest point of \code{X};
    \item a nearest neighbour map giving, for each pixel,
    the identifier of the nearest point of \code{X}.
  }
  If the argument \code{k} is specified, then the \code{k}-th nearest
  neighbours will be found.

  If \code{what="dist"} then only the distance map is returned.
  If \code{what="which"} then only the nearest neighbour map
  is returned.

  The argument \code{k} may be an integer or an integer vector.
  If it is a single integer, then the \code{k}-th nearest neighbours
  are computed. If it is a vector, then the \code{k[i]}-th nearest
  neighbours are computed for each entry \code{k[i]}. For example, setting
  \code{k=1:3} will compute the nearest, second-nearest and
  third-nearest neighbours. 
}
\section{Sorting data and pre-sorted data}{
  Read this section if you care about the speed of computation.
  
  For efficiency, the algorithm sorts the point pattern \code{X}
  into increasing order of the \eqn{x} coordinate
  or increasing order of the the \eqn{y} coordinate.
  Sorting is only an intermediate step;
  it does not affect the output, which is always given in the same
  order as the original data.
  
  By default (if \code{sortby="range"}),
  the sorting will occur on the coordinate that has the larger range of
  values (according to the frame of the enclosing window of \code{X}).
  If \code{sortby = "var"}), sorting will occur on the coordinate that
  has the greater variance (in the pattern \code{X}).
  Setting \code{sortby="x"} or \code{sortby = "y"} will specify that
  sorting should occur on the \eqn{x} or \eqn{y} coordinate, respectively.

  If the point pattern \code{X} is already
  sorted, then the argument \code{is.sorted.X}
  should be set to \code{TRUE}, and \code{sortby} should be set
  equal to \code{"x"} or \code{"y"} to indicate which coordinate
  is sorted.
}
\section{Warning About Ties}{
  Ties are possible: there may be two data points
  which lie exactly the same distance away from a particular pixel.
  This affects the results from \code{nnmap(what="which")}.
  The handling of ties is not well-defined: it is not consistent
  between different computers and different installations of \R.
  If there are ties, then different calls to \code{nnmap(what="which")}
  may give inconsistent results. For example, you may get a different answer
  from \code{nnmap(what="which",k=1)}
  and \code{nnmap(what="which", k=1:2)[[1]]}.
}
\value{
  A pixel image, or a list of pixel images.

  By default (if \code{what=c("dist", "which")}), the result is
  a list with two components \code{dist} and \code{which}
  containing the distance map and the nearest neighbour map.

  If \code{what="dist"} then the result is a real-valued pixel image
  containing the distance map.
  
  If \code{what="which"} then the result is an integer-valued pixel image
  containing the nearest neighbour map.

  If \code{k} is a vector of several integers, then the result is
  similar except that each pixel image is replaced by a list of
  pixel images, one for each entry of \code{k}.
}
\seealso{
  \code{\link{distmap}}
}
\examples{
  plot(nnmap(cells, 2, what="which"))
}
\author{
  \adrian
  
  ,
  \rolf
  ,
  and Jens Oehlschlaegel
}
\keyword{spatial}
\keyword{math}
