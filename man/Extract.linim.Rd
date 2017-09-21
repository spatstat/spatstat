\name{Extract.linim}
\alias{[.linim}
\title{Extract Subset of Pixel Image on Linear Network}
\description{
  Extract a subset of a pixel image on a linear network.
}
\usage{
  \method{[}{linim}(x, i, \dots, drop=TRUE)
}
\arguments{
  \item{x}{
    A pixel image on a linear network (object of class \code{"linim"}).
  }
  \item{i}{
    Spatial window defining the subregion.
    Either a spatial window (an object of class \code{"owin"}), or a
    logical-valued pixel image, 
    or any type of index that applies to a matrix,
    or a point pattern (an object of class \code{"lpp"} or \code{"ppp"}),
    or something that can be converted to a point pattern
    by \code{\link{as.lpp}} (using the network on which \code{x} is defined).
  }
  \item{\dots}{Additional arguments passed to \code{[.im}.}
  \item{drop}{Logical value indicating whether \code{NA} values should
  be omitted from the result.}
}
\value{
  Another pixel image on a linear network (object of class \code{"linim"})
  or a vector of pixel values.
}
\details{
  This function is a method for the subset operator \code{"["} for
  pixel images on linear networks (objects of class \code{"linim"}). 

  The pixel image \code{x} will be restricted to the 
  domain specified by \code{i}.

  Pixels outside the domain of \code{x} are assigned the value \code{NA};
  if \code{drop=TRUE} (the default) such \code{NA} values are deleted
  from the result; if \code{drop=FALSE}, then \code{NA} values are retained.

  If \code{i} is a window (or a logical-valued pixel image)
  then \code{x[i]} is another pixel image of class \code{"linim"},
  representing the restriction of \code{x} to the spatial domain
  specified by \code{i}.

  If \code{i} is a point pattern, then \code{x[i]} is the vector of
  pixel values of \code{x} at the locations specified by \code{i}.
}
\examples{
  M <- as.mask.psp(as.psp(simplenet))
  Z <- as.im(function(x,y){x}, W=M)
  Y <- linim(simplenet, Z)
  X <- runiflpp(4, simplenet)
  Y[X]
  Y[square(c(0.3, 0.6))]
}
\author{
  \adrian
}
\keyword{spatial}
\keyword{manip}
