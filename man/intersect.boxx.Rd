\name{intersect.boxx}
\alias{intersect.boxx}
\title{Intersection Of Boxes Of Arbitrary Dimension}
\description{
  Yields the intersection of boxes of arbitrary dimension (of class \code{"boxx"}).
}
\usage{
   intersect.boxx(\dots, fatal=FALSE)
}
\arguments{
  \item{\dots}{
    Boxes (of class \code{"boxx"}).
  }
  \item{fatal}{Logical.
    Determines what happens if the intersection is empty: If true
  }
}
\value{
  A box (object of class \code{"boxx"})
  or possibly \code{NULL}.
}
\details{
  If the intersection is empty, then if \code{fatal=FALSE}
  the result is \code{NULL},
  while if \code{fatal=TRUE} an error occurs.
}

\author{
\spatstatAuthors.
}
\seealso{
  \code{\link{intersect.owin}},
  \code{\link{boxx}}
}
\examples{
intersect.boxx(boxx(c(-1,1),c(0,2)), boxx(c(0,3),c(0,1)))
}
\keyword{spatial}
\keyword{math}
