\name{uniquemap.default}
\alias{uniquemap.default}
\alias{uniquemap.data.frame}
\alias{uniquemap.matrix}
\title{
  Map Duplicate Entries to Unique Entries
}
\description{
  Determine whether entries in a vector
  (or rows in a matrix or data frame)
  are duplicated,
  choose a unique representative for each set of duplicates,
  and map the duplicates to the unique representative.
}
\usage{
\method{uniquemap}{default}(x)

\method{uniquemap}{data.frame}(x)

\method{uniquemap}{matrix}(x)
}
\arguments{
  \item{x}{
    A vector, data frame or matrix,
    or another type of data.
  }
}
\details{
  The function \code{\link{uniquemap}} is generic, with methods
  for point patterns, data frames, and a default method.

  The default method expects a vector. It determines whether any entries
  of the vector \code{x} are duplicated,
  and constructs a mapping of the indices of \code{x}
  so that all duplicates are mapped to a unique representative index.

  The result is an integer vector \code{u} such that
  \code{u[j] = i} if
  the entries \code{x[i]} and \code{x[j]} are identical and
  point \code{i} has been chosen as the unique representative.
  The entry \code{u[i] = i} means either that point \code{i} is
  unique, or that it has been chosen as the unique representative
  of its equivalence class.

  The method for \code{data.frame} determines whether any rows of the
  data frame \code{x} are duplicated, and constructs a mapping of the
  row indices so that all duplicate rows are mapped to a unique
  representative row.
}
\value{
  An integer vector.
}
\author{
  \spatstatAuthors.
}
\seealso{
  \code{\link{uniquemap}},
  \code{\link{uniquemap.ppp}}
}
\examples{
  x <- c(3, 5, 2, 4, 2, 3)
  uniquemap(x)

  df <- data.frame(A=x, B=42)
  uniquemap(df)

  z <- cbind(x, 10-x)
  uniquemap(z)
}
\keyword{spatial}
\keyword{methods}
