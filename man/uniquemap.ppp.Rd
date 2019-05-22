\name{uniquemap.ppp}
\alias{uniquemap}
\alias{uniquemap.ppp}
\alias{uniquemap.lpp}
\alias{uniquemap.ppx}
\title{
  Map Duplicate Entries to Unique Entries
}
\description{
  Determine whether points in a point pattern are duplicated,
  choose a unique representative for each set of duplicates,
  and map the duplicates to the unique representative.
}
\usage{
uniquemap(x)

\method{uniquemap}{ppp}(x)

\method{uniquemap}{lpp}(x)

\method{uniquemap}{ppx}(x)
}
\arguments{
  \item{x}{
    A point pattern (object of class \code{"ppp"}, \code{"lpp"},
    \code{"pp3"} or \code{"ppx"}).
  }
}
\details{
  The function \code{uniquemap} is generic, with methods for
  point patterns and data frames.

  This function determines whether any points of \code{x} are
  duplicated, and constructs a mapping of the indices of \code{x}
  so that all duplicates are mapped to a unique representative index.

  The result is an integer vector \code{u} such that
  \code{u[j] = i} if
  the points \code{x[i]} and \code{x[j]} are identical and
  point \code{i} has been chosen as the unique representative.
  The entry \code{u[i] = i} means either that point \code{i} is
  unique, or that it has been chosen as the unique representative
  of its equivalence class.
}
\value{
  An integer vector.
}
\author{
  \spatstatAuthors.
}
\seealso{
  \code{\link{unique.ppp}}, \code{\link{duplicated.ppp}},
  \code{\link{uniquemap.default}}
}
\examples{
  Y <- runifpoint(4)
  X <- Y[c(1,2,3,4,2,1)]
  uniquemap(X)
}
\keyword{spatial}
\keyword{methods}
