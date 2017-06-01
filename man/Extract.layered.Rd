\name{Extract.layered}
\alias{[.layered}
\alias{[<-.layered}
\alias{[[<-.layered}
\title{Extract or Replace Subset of a Layered Object}
\description{
  Extract or replace some or all of the layers of a layered object,
  or extract a spatial subset of each layer.
}
\usage{
  \method{[}{layered}(x, i, j, drop=FALSE, ...)

  \method{[}{layered}(x, i) <- value

  \method{[[}{layered}(x, i) <- value
}
\arguments{
  \item{x}{
    A layered object (class \code{"layered"}).
  }
  \item{i}{
    Subset index for the list of layers.
    A logical vector, integer vector or character vector
    specifying which layers are to be extracted or replaced.
  }
  \item{j}{
    Subset index to be applied to the data in each layer.
    Typically a spatial window (class \code{"owin"}).
  }
  \item{drop}{
    Logical. If \code{i} specifies only a single layer
    and \code{drop=TRUE}, then the contents of this layer
    will be returned.
  }
  \item{\dots}{
     Additional arguments, passed to other subset methods
     if the subset index is a window.
  }
  \item{value}{List of objects which shall replace the designated
    subset, or an object which shall replace the designated element.
  }
}
\value{
  Usually an object of class \code{"layered"}.
}
\details{
  A layered object represents data that should be plotted in
  successive layers, for example, a background and a foreground.
  See \code{\link{layered}}.

  The function \code{[.layered}
  extracts a designated subset of a layered object.
  It is a method for \code{\link{[}} for the
  class \code{"layered"}.

  The functions \code{[<-.layered} and \code{[[<-.layered}
  replace a designated subset or designated entry of the object by new
  values. They are methods for \code{\link{[<-}} and \code{\link{[[<-}}
  for the \code{"layered"} class.

  The index \code{i} specifies which layers will be retained.
  It should be a valid subset index for the list of layers.

  The index \code{j} will be applied to each layer. It is typically
  a spatial window (class \code{"owin"}) so that each of the layers
  will be restricted to the same spatial region.
  Alternatively \code{j} may be any subset index
  which is permissible for the \code{"["} method for each of the layers.
}
\seealso{
  \code{\link{layered}}
}
\examples{
 D <- distmap(cells)
 L <- layered(D, cells,
              plotargs=list(list(ribbon=FALSE), list(pch=16)))

 L[-2]
 L[, square(0.5)]

 L[[3]] <- japanesepines
 L
}
\author{
  \spatstatAuthors
}
\keyword{spatial}
\keyword{manip}
