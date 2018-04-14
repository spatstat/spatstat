\name{tess}
\alias{tess}
\title{Create a Tessellation}
\description{
  Creates an object of class \code{"tess"} representing a tessellation
  of a spatial region.
}
\usage{
  tess(..., xgrid = NULL, ygrid = NULL, tiles = NULL, image = NULL,
            window=NULL, marks=NULL, keepempty=FALSE, unitname=NULL, check=TRUE)
}
\arguments{
  \item{\dots}{Ignored.}
  \item{xgrid,ygrid}{Cartesian coordinates of vertical and
    horizontal lines determining a grid of rectangles.
    Incompatible with other arguments.
  }
  \item{tiles}{List of tiles in the tessellation.
    A list, each of whose elements is a window
    (object of class \code{"owin"}). Incompatible with other arguments.
  }
  \item{image}{
    Pixel image which specifies the tessellation.
    Incompatible with other arguments.
  }
  \item{window}{
    Optional.
    The spatial region which is tessellated (i.e. the union of all the tiles).
    An object of class \code{"owin"}.
  }
  \item{marks}{
    Optional vector or data frame of marks associated with the tiles.
  }
  \item{keepempty}{
    Logical flag indicating whether empty tiles should be retained
    or deleted.
  }
  \item{unitname}{
    Optional. Name of unit of length. Either a single character string,
    or a vector of two character strings giving the
    singular and plural forms, respectively.
    If this argument is missing or \code{NULL},
    information about the unitname will be
    extracted from the other arguments.
    If this argument is given, it overrides any other information
    about the unitname.
  }
  \item{check}{
    Logical value indicating whether to check the validity of the
    input data. It is strongly recommended to use the default
    value \code{check=TRUE}.
  }
}
\details{
  A tessellation is a collection of disjoint spatial regions
  (called \emph{tiles}) that fit together to form a larger spatial
  region. This command creates an object of class \code{"tess"} that
  represents a tessellation.

  Three types of tessellation are supported:
  \describe{
    \item{rectangular:}{
      tiles are rectangles, with sides parallel to the \code{x} and
      \code{y} axes. They may or may not have equal size and shape.
      The arguments \code{xgrid} and \code{ygrid} determine
      the positions of the vertical and horizontal grid lines,
      respectively. (See \code{\link{quadrats}} for another way to do this.)
    }
    \item{tile list:}{
      tiles are arbitrary spatial regions.
      The argument \code{tiles} is a list of these tiles,
      which are objects of class \code{"owin"}.
    }
    \item{pixel image:}{
      Tiles are subsets of a fine grid of pixels.
      The argument \code{image} is a pixel image (object of class
      \code{"im"}) with factor values. Each level of the factor
      represents a different tile of the tessellation. The pixels that
      have a particular value of the factor constitute a tile. 
    }
  }

  The optional argument \code{window} specifies the spatial region
  formed by the union of all the tiles. In other words it specifies the
  spatial region that is divided into tiles by the tessellation.
  If this argument is missing or \code{NULL}, it will be determined by
  computing the set union of all the tiles. This is a time-consuming
  computation. For efficiency it is advisable to specify the window.
  Note that the validity of the window will not be checked.

  Empty tiles may occur, either because one of the entries in the list
  \code{tiles} is an empty window, or because one of the levels of the
  factor-valued pixel image \code{image} does not occur in the pixel data.
  When \code{keepempty=TRUE}, empty tiles are permitted. 
  When \code{keepempty=FALSE} (the default), tiles are not allowed to be
  empty, and any empty tiles will be removed from the tessellation.

  There are methods for \code{print}, \code{plot}, \code{[} and \code{[<-}
  for tessellations. Use \code{\link{tiles}} to extract the list of
  tiles in a tessellation, \code{\link{tilenames}} to extract the names
  of the tiles, and \code{\link{tile.areas}} to compute their
  areas.

  The tiles may have marks, which can be extracted by
  \code{\link{marks.tess}} and changed by \code{\link{marks<-.tess}}.

  Tessellations can be used to classify the points of
  a point pattern, in \code{\link{split.ppp}}, \code{\link{cut.ppp}} and
  \code{\link{by.ppp}}.

  To construct particular tessellations, see
  \code{\link{quadrats}}, \code{\link{hextess}}, 
  \code{\link{dirichlet}}, \code{\link{delaunay}}
  and \code{\link{rpoislinetess}}.
}
\value{
  An object of class \code{"tess"} representing the tessellation.
}
\seealso{
  \code{\link{marks.tess}},
  \code{\link{plot.tess}},
  \code{\link{[.tess}},
  \code{\link{as.tess}},
  \code{\link{tiles}},
  \code{\link{intersect.tess}},
  \code{\link{split.ppp}},
  \code{\link{cut.ppp}},
  \code{\link{by.ppp}},
  \code{\link{bdist.tiles}},
  \code{\link{tile.areas}}.

  To construct particular tessellations, see
  \code{\link{quadrats}}, \code{\link{hextess}},
  \code{\link{dirichlet}}, \code{\link{delaunay}}
  and \code{\link{rpoislinetess}}.

  To divide space into pieces containing equal
  amounts of stuff, use \code{\link{quantess}}.
}
\examples{
  A <- tess(xgrid=0:4,ygrid=0:4)
  A
  plot(A)
  B <- A[c(1, 2, 5, 7, 9)]
  B
  v <- as.im(function(x,y){factor(round(5 * (x^2 + y^2)))}, W=owin())
  levels(v) <- letters[seq(length(levels(v)))]
  E <- tess(image=v)
  plot(E)
  G <- tess(image=v, marks=toupper(levels(v)), unitname="km")
  G
}
\author{
  \spatstatAuthors.
}
\keyword{spatial}
\keyword{datagen}
