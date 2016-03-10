\name{ppx}
\Rdversion{1.1}
\alias{ppx}
\title{
  Multidimensional Space-Time Point Pattern
}
\description{
  Creates a multidimensional space-time point pattern
  with any kind of coordinates and marks.
}
\usage{
ppx(data, domain=NULL, coord.type=NULL, simplify=FALSE)
}
\arguments{
  \item{data}{
    The coordinates and marks of the points.
    A \code{data.frame} or \code{hyperframe}.
  }
  \item{domain}{
    Optional. The space-time domain containing the points.
    An object in some appropriate format, or \code{NULL}.
  }
  \item{coord.type}{
    Character vector specifying how each column of
    \code{data} should be interpreted: as a spatial coordinate,
    a temporal coordinate, a local coordinate or a mark.
    Entries are partially matched to the values
    \code{"spatial"},
    \code{"temporal"},
    \code{"local"} and
    \code{"mark"}.
  }
  \item{simplify}{
    Logical value indicating whether to simplify the result
    in special cases.
    If \code{simplify=TRUE}, a two-dimensional point pattern
    will be returned as an object of class \code{"ppp"},
    and a three-dimensional point pattern will be returned as
    an object of class \code{"pp3"}.
    If \code{simplify=FALSE} (the default) then the result is always
    an object of class \code{"ppx"}.
  }
}
\details{
  An object of class \code{"ppx"}
  represents a marked point pattern
  in multidimensional space and/or time. There may be any
  number of spatial coordinates, any number of temporal coordinates,
  any number of local coordinates,
  and any number of mark variables. The individual marks may be
  atomic (numeric values, factor values, etc) or objects of any kind.

  The argument \code{data} should contain the coordinates and marks of
  the points. It should be a \code{data.frame} or more generally a
  \code{hyperframe} (see \code{\link{hyperframe}}) with
  one row of data for each point.

  Each column of \code{data} is either
  a spatial coordinate, a temporal coordinate,
  a local coordinate, or a mark variable.
  The argument \code{coord.type} determines how each column is interpreted.
  It should be a character vector, of length equal to the number of
  columns of \code{data}. It should contain strings that partially match
  the values \code{"spatial"}, \code{"temporal"}, \code{"local"} and
  \code{"mark"}. (The first letters will be sufficient.)
  
  By default (if \code{coord.type} is missing or \code{NULL}),
  columns of numerical data are assumed to represent
  spatial coordinates, while other columns are assumed to be marks.
}
\value{
  Usually an object of class \code{"ppx"}.
  If \code{simplify=TRUE} the result may be an object of
  class \code{"ppp"} or \code{"pp3"}.
}
\author{Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  and Rolf Turner \email{r.turner@auckland.ac.nz}
  
}
\seealso{
  \code{\link{pp3}},
  \code{\link{print.ppx}} 
}
\examples{
   df <- data.frame(x=runif(4),y=runif(4),t=runif(4),
                    age=rep(c("old", "new"), 2),
                    size=runif(4))
   X <- ppx(data=df, coord.type=c("s","s","t","m","m"))
   X

   val <- 20 * runif(4)
   E <- lapply(val, function(s) { rpoispp(s) })
   hf <- hyperframe(t=val, e=as.listof(E))
   Z <- ppx(data=hf, domain=c(0,1))
   Z
}
\keyword{spatial}
\keyword{datagen}
