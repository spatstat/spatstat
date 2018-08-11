\name{Extract.hyperframe}
\alias{[.hyperframe}
\alias{[<-.hyperframe}
\alias{$.hyperframe}
\alias{$<-.hyperframe}
\title{Extract or Replace Subset of Hyperframe}
\description{
  Extract or replace a subset of a hyperframe.
}
\usage{
  \method{[}{hyperframe}(x, i, j, drop, strip=drop, \dots)
  \method{[}{hyperframe}(x, i, j) <- value
  \method{$}{hyperframe}(x, name)
  \method{$}{hyperframe}(x, name) <- value
}
\arguments{
  \item{x}{
    A hyperframe (object of class \code{"hyperframe"}).
  }
  \item{i,j}{
    Row and column indices.
  }
  \item{drop,strip}{
     Logical values indicating what to do when the
     hyperframe has only one row or column.
     See Details.
  }
  \item{\dots}{
    Ignored. 
  }
  \item{name}{
    Name of a column of the hyperframe.
  }
  \item{value}{
    Replacement value for the subset. A hyperframe
    or (if the subset is a single column) a list or an atomic vector.
  }
}
\value{
  A hyperframe (of class \code{"hyperframe"}).
}
\details{
  These functions extract a designated subset of a hyperframe,
  or replace the designated subset with another hyperframe.

  The function \code{[.hyperframe} is a method for the subset operator
  \code{\link{[}} for the
  class \code{"hyperframe"}. It extracts the subset of \code{x}
  specified by the row index \code{i} and column index \code{j}.

  The argument \code{drop} determines whether the array structure
  will be discarded if possible. The argument \code{strip} determines
  whether the list structure in a row or column or cell will be discarded if
  possible.
  If \code{drop=FALSE} (the default), the return value is always
  a hyperframe or data frame.
  If \code{drop=TRUE}, and if the selected subset has only one row,
  or only one column, or both, then
  \itemize{
    \item{
       if \code{strip=FALSE}, the result is a
       list, with one entry for each array cell that was selected.
    }
    \item{
       if \code{strip=TRUE},
       \itemize{
         \item if the subset has one row containing several columns,
         the result is a list or (if possible) an atomic vector;
         \item if the subset has one column containing several rows,
         the result is a list or (if possible) an atomic vector;
         \item if the subset has exactly one row and exactly one column,
         the result is the object (or atomic value)
         contained in this row and column.
       }
    }
  }
  
  The function \code{[<-.hyperframe} is a method for the
  subset replacement operator \code{\link{[<-}} for the
  class \code{"hyperframe"}. It replaces the designated
  subset with the hyperframe \code{value}.
  The subset of \code{x} to be replaced is designated by
  the arguments \code{i} and \code{j} as above.
  The replacement \code{value} should be a hyperframe with the
  appropriate dimensions, or (if the specified subset is a single
  column) a list of the appropriate length.

  The function \code{$.hyperframe} is a method for \code{\link{$}}
  for hyperframes. It extracts the relevant column of the hyperframe.
  The result is always a list (i.e. equivalent to using
  \code{[.hyperframe} with \code{strip=FALSE}).

  The function \code{$<-.hyperframe} is a method for \code{\link{$<-}}
  for hyperframes. It replaces the relevant column of the hyperframe.
  The replacement value should be a list of the appropriate length.
}
\seealso{
  \code{\link{hyperframe}}
}
\examples{
  h <- hyperframe(X=list(square(1), square(2)), Y=list(sin, cos))
  h
  h[1, ]
  h[1, ,drop=TRUE]
  h[ , 1]
  h[ , 1, drop=TRUE]
  h[1,1]
  h[1,1,drop=TRUE]
  h[1,1,drop=TRUE,strip=FALSE]
  h[1,1] <- list(square(3))
  # extract column
  h$X
  # replace existing column
  h$Y <- list(cells, cells)
  # add new column
  h$Z <- list(cells, cells)
}
\author{
  \spatstatAuthors
}
\keyword{spatial}
\keyword{manip}
