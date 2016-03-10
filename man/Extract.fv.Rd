\name{Extract.fv}
\alias{[.fv}
\alias{[<-.fv}
\alias{$<-.fv}
\title{Extract or Replace Subset of Function Values}
\description{
  Extract or replace a subset of an object of class \code{"fv"}.
}
\usage{
  \method{[}{fv}(x, i, j, \dots, drop=FALSE)

  \method{[}{fv}(x, i, j) <- value

  \method{$}{fv}(x, name) <- value
}
\arguments{
  \item{x}{
    a function value object,
    of class \code{"fv"} (see \code{\link{fv.object}}).
    Essentially a data frame.
  }
  \item{i}{
    any appropriate subset index.
    Selects a subset of the rows of the data frame, i.e.
    a subset of the domain of the function(s) represented by \code{x}.
  }
  \item{j}{
    any appropriate subset index for the columns of the data frame.
    Selects some of the functions present in \code{x}.
  }
  \item{name}{
    the name of a column of the data frame.
  }
  \item{\dots}{
    Ignored.
  }
  \item{drop}{
     Logical. If \code{TRUE}, the result is a data frame or vector
     containing the selected rows and columns of data.
     If \code{FALSE} (the default),
     the result is another object of class \code{"fv"}.
  }
  \item{value}{
     Replacement value for the column or columns selected by
     \code{name} or \code{j}.
  }
} 
\value{
  The result of \code{[.fv} with \code{drop=TRUE}
  is a data frame or vector. 

  Otherwise, the result is another object of class \code{"fv"}.
}
\details{
  These functions extract a designated subset of an object of class
  \code{"fv"}, or replace the designated subset with other data,
  or delete the designated subset.

  The subset is specified by the 
  row index \code{i} and column index \code{j}, or
  by the column name \code{name}. Either \code{i} or \code{j}
  may be missing, or both may be missing.

  The function \code{[.fv} is a method for the generic operator
  \code{\link{[}} for the class \code{"fv"}. It extracts the
  designated subset of \code{x}, and returns it as
  another object of class \code{"fv"} (if \code{drop=FALSE})
  or as a data frame or vector (if \code{drop=TRUE}).

  The function \code{[<-.fv} is a method for the generic operator
  \code{\link{[<-}} for the class \code{"fv"}.
  If \code{value} is \code{NULL}, the designated subset of \code{x} will be
  deleted from \code{x}.
  Otherwise, the designated subset of \code{x} will be
  replaced by the data contained in \code{value}.
  The return value is the modified object \code{x}.

  The function \code{$<-.fv} is a method for the generic operator
  \code{\link{$<-}} for the class \code{"fv"}.
  If \code{value} is \code{NULL}, the designated column of \code{x} will be
  deleted from \code{x}.
  Otherwise, the designated column of \code{x} will be
  replaced by the data contained in \code{value}.
  The return value is the modified object \code{x}.
}
\seealso{
  \code{\link{fv.object}}
}
\examples{
 K <- Kest(cells)

 # discard the estimates of K(r) for r  > 0.1
 Ksub <- K[K$r <= 0.1, ]

 # extract the border method estimates
 bor <- K[ , "border", drop=TRUE]
 # or equivalently
 bor <- K$border

 # remove the border-method estimates
 K$border <- NULL
 K
}
\author{
  Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  and
  Rolf Turner \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{manip}
