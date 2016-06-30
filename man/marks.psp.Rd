\name{marks.psp}
\alias{marks.psp}
\alias{marks<-.psp}
\title{Marks of a Line Segment Pattern}
\description{
  Extract or change the marks attached to
  a line segment pattern.
}
\usage{
\method{marks}{psp}(x, \dots, dfok=TRUE)           
\method{marks}{psp}(x, \dots) <- value
}
\arguments{
  \item{x}{
    Line segment pattern dataset (object of class \code{"psp"}).
  }
  \item{\dots}{
    Ignored.
  }
  \item{dfok}{
    Logical. If \code{FALSE}, data frames of marks are not permitted
    and will generate an error.
  }
  \item{value}{
    Vector or data frame of mark values,
    or \code{NULL}.
  }
}
\value{
  For \code{marks(x)}, the result is a vector, factor or data frame,
  containing the mark values attached to the line segments of \code{x}.
  If there are no marks, the result is \code{NULL}.

  For \code{marks(x) <- value}, the result is the updated line segment pattern
  \code{x} (with the side-effect that the dataset \code{x} is updated in
  the current environment).
}
\details{
  These functions extract or change the marks
  attached to each of the line segments in the pattern \code{x}.
  They are methods for the generic functions
  \code{\link{marks}} and \code{\link{marks<-}}
  for the class \code{"psp"} of line segment patterns.
    
  The expression \code{marks(x)} extracts the marks of \code{x}.
  The assignment \code{marks(x) <- value} assigns new marks to the
  dataset \code{x}, and updates the dataset \code{x} in the current
  environment. 
  
  The marks can be a vector, a factor, or a data frame.
  
  For the assignment \code{marks(x) <- value}, the \code{value}
  should be a vector or factor of length equal to the number of
  segments in \code{x}, or a data frame with as many rows
  as there are segments in \code{x}. If \code{value} is a single value,
  or a data frame with one row, then it will be replicated
  so that the same marks will be attached to each segment.
  
  To remove marks, use \code{marks(x) <- NULL} or \code{unmark(x)}.
}
\seealso{
  \code{\link{psp.object}},
  \code{\link{marks}},
  \code{\link{marks<-}}
}
\examples{
  m <- data.frame(A=1:10, B=letters[1:10])
  X <- psp(runif(10), runif(10), runif(10), runif(10), window=owin(), marks=m)

  marks(X) 
  marks(X)[,2]
  marks(X) <- 42
  marks(X) <- NULL
}
\author{
  \spatstatAuthors.
}
\keyword{spatial}
\keyword{manip}
