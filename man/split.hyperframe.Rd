\name{split.hyperframe}
\alias{split.hyperframe}
\alias{split<-.hyperframe}
\title{
  Divide Hyperframe Into Subsets and Reassemble
}
\description{
  \code{split} divides the data \code{x} into subsets defined
  by \code{f}. The replacement form replaces values corresponding to
  such a division. 
}
\usage{
  \method{split}{hyperframe}(x, f, drop = FALSE, ...)

  \method{split}{hyperframe}(x, f, drop = FALSE, ...) <- value
}
\arguments{
  \item{x}{
    Hyperframe (object of class \code{"hyperframe"}).
  }
  \item{f}{
    a \code{factor} in the sense that \code{as.factor(f)} defines the
    grouping, or a list of such factors in which case their
    interaction is used for the grouping.
  }
  \item{drop}{
    logical value, indicating whether levels that do not occur should be
    dropped from the result.
  }
  \item{value}{
    a list of hyperframes which arose (or could have arisen)
    from the command \code{split(x,f,drop=drop)}.
  }
  \item{\dots}{
    Ignored.
  }
}
\details{
  These are methods for the generic functions
  \code{\link{split}} and \code{\link{split<-}}
  for hyperframes (objects of class \code{"hyperframe"}).

  A hyperframe is like a data frame, except that its entries
  can be objects of any kind. The behaviour of these methods
  is analogous to the corresponding methods for data frames.
}
\value{
  The value returned from \code{split.hyperframe} is a list of
  hyperframe containing
  the values for the groups.  The components of the list are named
  by the levels of \code{f} (after converting to a factor, or if already
  a factor and \code{drop = TRUE}, dropping unused levels).

  The replacement method \code{split<-.hyperframe} returns
  a new hyperframe \code{x} for which \code{split(x,f)} equals \code{value}.
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/},
  Rolf Turner
  \email{r.turner@auckland.ac.nz}
  and Ege Rubak
}
\seealso{
  \code{\link{hyperframe}}, \code{\link{[.hyperframe}}
}
\examples{
   split(pyramidal, pyramidal$group)
}
\keyword{spatial}
\keyword{manip}
