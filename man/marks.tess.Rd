\name{marks.tess}
\alias{marks.tess}
\alias{marks<-.tess}
\alias{unmark.tess}
\title{Marks of a Tessellation}
\description{
  Extract or change the marks attached to
  the tiles of a tessellation.
}
\usage{
\method{marks}{tess}(x, \dots)           

\method{marks}{tess}(x, \dots) <- value

\method{unmark}{tess}(X)
}
\arguments{
  \item{x,X}{
    Tessellation (object of class \code{"tess"}).
  }
  \item{\dots}{
    Ignored.
  }
  \item{value}{
    Vector or data frame of mark values,
    or \code{NULL}.
  }
}
\value{
  For \code{marks(x)}, the result is a vector, factor or data frame,
  containing the mark values attached to the tiles of \code{x}.
  If there are no marks, the result is \code{NULL}.

  For \code{unmark(x)}, the result is the tessellation without marks.
  
  For \code{marks(x) <- value}, the result is the updated tessellation
  \code{x} (with the side-effect that the dataset \code{x} is updated in
  the current environment).
}
\details{
  These functions extract or change the marks
  attached to each of the tiles in the tessellation \code{x}.
  They are methods for the generic functions
  \code{\link{marks}} and \code{\link{marks<-}}
  for the class \code{"tess"} of tessellations.
    
  The expression \code{marks(x)} extracts the marks of \code{x}.
  The assignment \code{marks(x) <- value} assigns new marks to the
  dataset \code{x}, and updates the dataset \code{x} in the current
  environment. 
  
  The marks can be a vector, a factor, or a data frame.
  
  For the assignment \code{marks(x) <- value}, the \code{value}
  should be a vector or factor of length equal to the number of
  tiles in \code{x}, or a data frame with as many rows
  as there are tiles in \code{x}. If \code{value} is a single value,
  or a data frame with one row, then it will be replicated
  so that the same marks will be attached to each tile.
  
  To remove marks, use \code{marks(x) <- NULL} or \code{unmark(x)}.
}
\seealso{
  \code{\link{marks}},
  \code{\link{marks<-}}
}
\examples{
  D <- dirichlet(cells)
  marks(D) <- tile.areas(D)
}
\author{
Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}


Rolf Turner \email{r.turner@auckland.ac.nz}

and Ege Rubak \email{rubak@math.aau.dk}

}
\keyword{spatial}
\keyword{manip}
