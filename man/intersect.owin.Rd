\name{intersect.owin}
\alias{intersect.owin}
\alias{union.owin}
\alias{setminus.owin}
\title{Intersection, Union or Set Subtraction of Two Windows}
\description{
  Yields the intersection, union or set subtraction of two windows.
}
\usage{
   intersect.owin(A, B, \dots, fatal=TRUE)
   union.owin(A,B, \dots)
   setminus.owin(A,B, \dots)
}
\arguments{
  \item{A}{A window object (see Details).}
  \item{B}{A window object.}
  \item{\dots}{
    Optional arguments passed to \code{\link{as.mask}}
    to control the discretisation, if required.
  }
  \item{fatal}{Logical.
    Determines what happens if the intersection is empty.
  }
}
\value{
  A window (object of class \code{"owin"}).
}
\details{
  The function \code{intersect.owin} computes the intersection between the
  two windows \code{A} and \code{B}, while
  \code{union.owin} computes their union.
  The function \code{setminus.owin} computes the intersection of
  \code{A} with the complement of \code{B}.

  The arguments \code{A} and \code{B} must be window objects
  (either objects of class \code{"owin"}, or data that can be
  coerced to this class by \code{\link{as.owin}}).

  If the intersection is empty, then if \code{fatal=FALSE}
  the result is NULL, while if \code{fatal=TRUE} an error occurs.

  The intersection or union of more than two windows can also be
  computed. For \code{intersect.owin} and \code{union.owin} the
  arguments \code{\dots} can include additional window objects.
}

\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\seealso{
  \code{\link{is.subset.owin}},
  \code{\link{overlap.owin}},
  \code{\link{bounding.box}},
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
   plot(bounding.box(A,B), main="intersection")
   w <- intersect.owin(A, B)
   plot(w, add=TRUE, col="lightblue")
   plot(A, add=TRUE)
   plot(B, add=TRUE)

   plot(bounding.box(A,B), main="union")
   w <- union.owin(A,B)
   plot(w, add=TRUE, col="lightblue")   
   plot(A, add=TRUE)
   plot(B, add=TRUE)

   plot(bounding.box(A,B), main="set minus")
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
