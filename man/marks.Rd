\name{marks}
\alias{marks}
\alias{marks.ppp}           
\alias{marks.ppx}           
\alias{marks<-}           
\alias{marks<-.ppp}           
\alias{marks<-.ppx}           
\alias{setmarks}
\alias{\%mark\%}  %DoNotExport 
%NAMESPACE export("%mark%")
\title{Marks of a Point Pattern}
\description{
  Extract or change the marks attached to
  a point pattern dataset.
}
\usage{
marks(x, \dots)

\method{marks}{ppp}(x, \dots, dfok=TRUE, drop=TRUE)           

\method{marks}{ppx}(x, \dots, drop=TRUE)

marks(x, \dots) <- value

\method{marks}{ppp}(x, \dots, dfok=TRUE, drop=TRUE) <- value

\method{marks}{ppx}(x, \dots) <- value

setmarks(x, value)

x \%mark\% value
}
\arguments{
  \item{x}{
    Point pattern dataset (object of class \code{"ppp"} or \code{"ppx"}).
  }
  \item{\dots}{
    Ignored.
  }
  \item{dfok}{
    Logical. If \code{FALSE}, data frames of marks are not permitted
    and will generate an error.
  }
  \item{drop}{
    Logical. If \code{TRUE}, a data frame consisting of a single column
    of marks will be converted to a vector or factor.
  }
  \item{value}{
    Replacement value. A vector, data frame or hyperframe of mark values,
    or \code{NULL}. 
  }
}
\value{
  For \code{marks(x)}, the result is a vector, factor, data frame or hyperframe,
  containing the mark values attached to the points of \code{x}.

  For \code{marks(x) <- value}, the result is the updated point pattern
  \code{x} (with the side-effect that the dataset \code{x} is updated in
  the current environment).

  For \code{setmarks(x,value)} and \code{x \%mark\% value}, the return value
  is the point pattern obtained by replacing the
  marks of \code{x} by \code{value}.
}
\details{
  These functions extract or change the marks
  attached to the points of the point pattern \code{x}.
  
  The expression \code{marks(x)} extracts the marks of \code{x}.
  The assignment \code{marks(x) <- value} assigns new marks to the
  dataset \code{x}, and updates the dataset \code{x} in the current
  environment. The expression \code{setmarks(x,value)}
  or equivalently \code{x \%mark\% value} returns a point pattern
  obtained by replacing the marks of \code{x} by \code{value}, but does
  not change the dataset \code{x} itself.
  
  For point patterns in two-dimensional space (objects of class
  \code{"ppp"}) the marks can be a vector, a factor, or a data frame.
  
  For general point patterns (objects of class "ppx") the
  marks can be a vector, a factor, a data frame or a
  hyperframe.
  
  For the assignment \code{marks(x) <- value}, the \code{value}
  should be a vector or factor of length equal to the number of
  points in \code{x}, or a data frame or hyperframe with as many rows
  as there are points in \code{x}. If \code{value} is a single value,
  or a data frame or hyperframe with one row, then it will be replicated
  so that the same marks will be attached to each point.
  
  To remove marks, use \code{marks(x) <- NULL} or
  \code{\link{unmark}(x)}.
  
  Use \code{\link{ppp}} or \code{\link{ppx}}
  to create point patterns in more general
  situations.
}
\seealso{
  \code{\link{ppp.object}},
  \code{\link{ppx}},
  \code{\link{unmark}},
  \code{\link{hyperframe}}
}
\examples{
   X <- amacrine
   # extract marks
   m <- marks(X)
   # recode the mark values "off", "on" as 0, 1
   marks(X) <- as.integer(m == "on")
}
\author{
  Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au},
 Rolf Turner \email{r.turner@auckland.ac.nz}
 and Ege Rubak \email{rubak@math.aau.dk}
}
\keyword{spatial}
\keyword{manip}
