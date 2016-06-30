\name{intersect.owin}
\alias{intersect.owin}
\alias{union.owin}
\alias{setminus.owin}
\title{Intersection, Union or Set Subtraction of Windows}
\description{
  Yields the intersection, union or set subtraction of windows.
}
\usage{
   intersect.owin(\dots, fatal=TRUE, p)
   union.owin(\dots, p)
   setminus.owin(A, B, \dots, p)
}
\arguments{
  \item{A,B}{Windows (objects of class \code{"owin"})
    or data acceptable to \code{\link{as.owin}}.}
  \item{\dots}{
    Windows,
    or arguments passed to \code{\link{as.mask}}
    to control the discretisation.
  }
  \item{fatal}{Logical.
    Determines what happens if the intersection is empty.
  }
  \item{p}{
    Optional list of parameters passed to
    \code{\link[polyclip]{polyclip}} to control the
    accuracy of polygon geometry.
  }
}
\value{
  A window (object of class \code{"owin"})
  or possibly \code{NULL}.
}
\details{
  The function \code{intersect.owin} computes the intersection between
  the windows given in \code{\dots}, while 
  \code{union.owin} computes their union.
  The function \code{setminus.owin} computes the intersection of
  \code{A} with the complement of \code{B}.

  For \code{intersect.owin} and \code{union.owin}, 
  the arguments \code{\dots} must be either
  \itemize{
    \item window objects of class \code{"owin"},
    \item data that can be coerced to this class by
    \code{\link{as.owin}}),
    \item lists of windows, of class \code{"solist"},
    \item named arguments of \code{\link{as.mask}} to control
    the discretisation if required.
  }
  For \code{setminus.owin}, the arguments \code{\dots}
  must be named arguments of \code{\link{as.mask}}.

  If the intersection is empty, then if \code{fatal=FALSE}
  the result is NULL, while if \code{fatal=TRUE} an error occurs.
}

\author{
\adrian


\rolf

and \ege

}
\seealso{
  \code{\link{is.subset.owin}},
  \code{\link{overlap.owin}},
  \code{\link{boundingbox}},
  \code{\link{owin.object}}
}
\examples{
# rectangles
   u <- unit.square()
   v <- owin(c(0.5,3.5), c(0.4,2.5))
# polygon
   data(letterR)
# mask
   m <- as.mask(letterR)

# two rectangles
   intersect.owin(u, v) 
   union.owin(u,v)
   setminus.owin(u,v)

# polygon and rectangle
   intersect.owin(letterR, v)
   union.owin(letterR,v)
   setminus.owin(letterR,v)

# mask and rectangle
   intersect.owin(m, v)
   union.owin(m,v)
   setminus.owin(m,v)

# mask and polygon
   p <- rotate(v, 0.2)
   intersect.owin(m, p)
   union.owin(m,p)
   setminus.owin(m,p)

# two polygons
   A <- letterR
   B <- rotate(letterR, 0.2)
   plot(boundingbox(A,B), main="intersection")
   w <- intersect.owin(A, B)
   plot(w, add=TRUE, col="lightblue")
   plot(A, add=TRUE)
   plot(B, add=TRUE)

   plot(boundingbox(A,B), main="union")
   w <- union.owin(A,B)
   plot(w, add=TRUE, col="lightblue")   
   plot(A, add=TRUE)
   plot(B, add=TRUE)

   plot(boundingbox(A,B), main="set minus")
   w <- setminus.owin(A,B)
   plot(w, add=TRUE, col="lightblue")   
   plot(A, add=TRUE)
   plot(B, add=TRUE)

# intersection and union of three windows
   C <- shift(B, c(0.2, 0.3))
   plot(union.owin(A,B,C))
   plot(intersect.owin(A,B,C))
}
\keyword{spatial}
\keyword{math}
