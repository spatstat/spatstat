\name{hyperframe}
\alias{hyperframe}
\title{Hyper Data Frame}
\description{
  Create a hyperframe: a two-dimensional array in which each column
  consists of values of the same atomic type (like the columns
  of a data frame) or objects of the same class. 
}
\usage{
   hyperframe(...,
             row.names=NULL, check.rows=FALSE, check.names=TRUE,
             stringsAsFactors=default.stringsAsFactors())
}
\arguments{
  \item{\dots}{
    Arguments of the form \code{value} or \code{tag=value}.
    Each \code{value} is either an atomic vector, or a list of objects of the
    same class, or a single atomic value, or a single object.
    Each \code{value} will become a column of the array.
    The \code{tag} determines the name of the column. See Details.
  }
  \item{row.names,check.rows,check.names,stringsAsFactors}{
    Arguments passed to \code{\link{data.frame}} controlling the
    names of the rows, whether to check that rows are consistent,
    whether to check validity of the column names, and whether to
    convert character columns to factors.
  }
}
\details{
  A hyperframe is like a data frame, except that its entries
  can be objects of any kind.

  A hyperframe is a two-dimensional array in which each column consists of
  values of one atomic type (as in a data frame) or consists of
  objects of one class.

  The arguments \code{\dots} are any number of arguments of
  the form \code{value} or \code{tag=value}. Each \code{value} will
  become a column of the array. The \code{tag} determines the name
  of the column.

  Each \code{value} can be either
  \itemize{
    \item an atomic vector or factor
    (i.e. numeric vector, integer vector, character vector, logical
    vector, complex vector or factor)
    \item a list of objects which are all of the same class
    \item one atomic value, which will be replicated to make an atomic
    vector or factor
    \item one object, which will be replicated to make a list of objects.
  }

  All columns (vectors, factors and lists) must be of the same length,
  if their length is greater than 1. 
}
\section{Methods for Hyperframes}{
  There are methods for
  \code{print}, \code{plot}, \code{summary}, \code{with}, \code{split},
  \code{[}, \code{[<},\code{$}, \code{$<-},
  \code{names}, \code{as.data.frame} \code{as.list},
  \code{cbind} and \code{rbind} for the class of hyperframes. There is also
  \code{is.hyperframe} and \code{\link{as.hyperframe}}.
}
\value{
  An object of class \code{"hyperframe"}.
}
\author{\adrian
  
  
  and \rolf
  
}
\seealso{
  \code{\link{as.hyperframe}},
  \code{\link{as.hyperframe.ppx}},
  \code{\link{plot.hyperframe}},
  \code{\link{[.hyperframe}},
  \code{\link{with.hyperframe}},
  \code{\link{split.hyperframe}},
  \code{\link{as.data.frame.hyperframe}},
  \code{\link{cbind.hyperframe}},
  \code{\link{rbind.hyperframe}}
}
\examples{
 # equivalent to a data frame
  hyperframe(X=1:10, Y=3)

 # list of functions
  hyperframe(f=list(sin, cos, tan))

 # table of functions and matching expressions
  hyperframe(f=list(sin, cos, tan),
             e=list(expression(sin(x)), expression(cos(x)), expression(tan(x))))

  hyperframe(X=1:10, Y=letters[1:10], Z=factor(letters[1:10]),
    stringsAsFactors=FALSE)

  lambda <- runif(4, min=50, max=100)
  X <- lapply(as.list(lambda), function(x) { rpoispp(x) })
  h <- hyperframe(lambda=lambda, X=X)
  h

  h$lambda2 <- lambda^2
  h[, "lambda3"] <- lambda^3
  h[, "Y"] <- X
}
\keyword{spatial}
\keyword{manip}
