\name{unstack.solist}
\alias{unstack.solist}
\alias{unstack.layered}
\title{
  Unstack Each Spatial Object in a List of Objects
}
\description{
  Given a list of two-dimensional spatial objects,
  apply 
  
}
\usage{
\method{unstack}{solist}(x, \dots)

\method{unstack}{layered}(x, \dots)
}
\arguments{
  \item{x}{
    An object of class \code{"solist"} or \code{"layered"}
    representing a list of two-dimensional spatial objects.
  }
  \item{\dots}{
    Ignored.
  }
}
\details{
  The functions defined here are methods for the generic
  \code{\link[utils]{unstack}}. They expect the argument \code{x}
  to be a list of spatial objects, of class \code{"solist"}
  or \code{"layered"}.

  Each spatial object in the list \code{x}
  will be unstacked by applying the relevant method for
  \code{\link[utils]{unstack}}. This means that
  \itemize{
    \item a marked point pattern with several columns of marks
    will be separated into several point patterns, each having
    a single column of marks
    \item a measure with \eqn{k}-dimensional vector values
    will be separated into \eqn{k} measures with scalar values
  }
  The resulting unstacked objects
  will be collected into a list of the same kind as \code{x}.
  Typically the length of \code{unstack(x)} is greater than
  the length of \code{x}.
}
\value{
  A list belonging to the same class as \code{x}.
}
\author{
  \spatstatAuthors.
}
\seealso{
  \code{\link[utils]{unstack}}

  \code{\link{unstack.msr}},
  \code{\link{unstack.ppp}},
  \code{\link{unstack.lpp}},
  \code{\link{unstack.psp}}
}
\examples{
  A <- solist(finpines=finpines, cells=cells)
  A
  unstack(A)
  B <- layered(fin=finpines, loc=unmark(finpines),
              plotargs=list(list(), list(pch=16)))
  B
  plot(B)
  unstack(B)
  plot(unstack(B))
}
\keyword{spatial}
\keyword{manip}
