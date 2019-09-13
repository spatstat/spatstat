\name{unstack.ppp}
\alias{unstack.ppp}
\alias{unstack.psp}
\alias{unstack.lpp}
\alias{unstack.tess}
\alias{unstack.lintess}
\title{
  Separate Multiple Columns of Marks 
}
\description{
  Given a spatial pattern with several columns of marks,
  take one column at a time, and return a list of spatial patterns
  each having only one column of marks.
}
\usage{
\method{unstack}{ppp}(x, \dots)

\method{unstack}{psp}(x, \dots)

\method{unstack}{lpp}(x, \dots)

\method{unstack}{tess}(x, \dots)

\method{unstack}{lintess}(x, \dots)
}
\arguments{
  \item{x}{
    A spatial point pattern
    (object of class \code{"ppp"} or \code{"lpp"})
    or a spatial pattern of line segments (object of class \code{"psp"})
    or a spatial tessellation (object of class \code{"tess"})
    or a tessellation on a linear network (object of class \code{"lintess"}).
  }
  \item{\dots}{
    Ignored.
  }
}
\details{
  The functions defined here are methods for the generic
  \code{\link[utils]{unstack}}. The functions expect a spatial object
  \code{x} which has several columns of marks; they separate the columns,
  and return a list of spatial objects, each having only one column of marks.
  
  If \code{x} has several columns of marks (i.e. \code{marks(x)} is a
  matrix, data frame or hyperframe with several columns),
  then \code{y <- unstack(x)} is a list of spatial objects, each of the same
  kind as \code{x}. The \code{j}th entry \code{y[[j]]} is equivalent to
  \code{x} except that it only includes
  the \code{j}th column of \code{marks(x)}.
  
  If \code{x} has no marks, or has only a single column of marks,
  the result is a list consisting of one entry, which is \code{x}.
}
\value{
  A list, of class \code{"solist"}, whose entries are objects of the
  same type as \code{x}.
}
\author{
  \spatstatAuthors.
}
\seealso{
  \code{\link[utils]{unstack}}

  \code{\link{unstack.msr}}

  See also methods for the generic \code{\link[base]{split}} such as
  \code{\link{split.ppp}}.
}
\examples{
   finpines
   unstack(finpines)
}
\keyword{spatial}
\keyword{manip}
