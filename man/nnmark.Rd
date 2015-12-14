\name{nnmark}
\alias{nnmark}
\title{
  Mark of Nearest Neighbour 
}
\description{
  Given a marked point pattern dataset \code{X}
  this function computes, for each desired location \code{y},
  the mark attached to the nearest neighbour of \code{y} in \code{X}.
  The desired locations \code{y} can be either a pixel grid
  or the point pattern \code{X} itself.
}
\usage{
nnmark(X, \dots, k = 1, at=c("pixels", "points"))
}
\arguments{
  \item{X}{
    A marked point pattern (object of class \code{"ppp"}).
  }
  \item{\dots}{
    Arguments passed to \code{\link{as.mask}} to determine the
    pixel resolution.
  }
  \item{k}{
    Single integer. The \code{k}th nearest data point will be used.
  }
  \item{at}{
    String specifying whether to compute the values
    at a grid of pixel locations (\code{at="pixels"}) or
    only at the points of \code{X} (\code{at="points"}).
  }
}
\details{
  Given a marked point pattern dataset \code{X}
  this function computes, for each desired location \code{y},
  the mark attached to the point of \code{X} that is nearest
  to \code{y}. The desired locations \code{y} can be either a pixel grid
  or the point pattern \code{X} itself.

  The argument \code{X} must be a marked point pattern (object
  of class \code{"ppp"}, see \code{\link{ppp.object}}).
  The marks are allowed to be a vector or a data frame.

  \itemize{
    \item
    If \code{at="points"}, then for each point in \code{X},
    the algorithm finds the nearest \emph{other} point in \code{X},
    and extracts the mark attached to it.
    The result is a vector or data frame containing the marks
    of the neighbours of each point.
    \item
    If \code{at="pixels"} (the default), then for each pixel
    in a rectangular grid, the algorithm finds the nearest point in \code{X},
    and extracts the mark attached to it.
    The result is an image or a list of images containing the marks
    of the neighbours of each pixel.
    The pixel resolution is controlled by the arguments \code{\dots}
    passed to \code{\link{as.mask}}.
  }

  If the argument \code{k} is given, then the \code{k}-th nearest
  neighbour will be used.
}
\value{
  \emph{If \code{X} has a single column of marks:}
  \itemize{
    \item 
    If \code{at="pixels"} (the default), the result is
    a pixel image (object of class \code{"im"}). 
    The value at each pixel is the mark attached
    to the nearest point of \code{X}.
    \item
    If \code{at="points"}, the result is a vector or factor
    of length equal to the number of points in \code{X}.
    Entries are the mark values of the
    nearest neighbours of each point of \code{X}.
  }
  \emph{If \code{X} has a data frame of marks:}
  \itemize{
    \item 
    If \code{at="pixels"} (the default), the result is a named list of 
    pixel images (object of class \code{"im"}). There is one
    image for each column of marks. This list also belongs to
    the class \code{"solist"}, for which there is a plot method.
    \item
    If \code{at="points"}, the result is a data frame
    with one row for each point of \code{X},
    Entries are the mark values of the
    nearest neighbours of each point of \code{X}.
  }
}
\author{Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  Rolf Turner \email{r.turner@auckland.ac.nz}
  
  and Ege Rubak \email{rubak@math.aau.dk}
  
}
\seealso{
  \code{\link{Smooth.ppp}},
  \code{\link{marktable}},
  \code{\link{nnwhich}}
}
\examples{
  plot(nnmark(ants))
  v <- nnmark(ants, at="points")
  v[1:10]
  plot(nnmark(finpines))
  vf <- nnmark(finpines, at="points")
  vf[1:5,]
}
\keyword{spatial}
\keyword{methods}
\keyword{smooth}
