\name{repairNetwork}
\alias{repairNetwork}
\title{
  Repair Internal Data in a Linear Network
}
\description{
  Detect and repair inconsistencies or duplication in the internal data
  of a network object.
}
\usage{
repairNetwork(X)
}
\arguments{
  \item{X}{
    A linear network (object of class \code{"linnet"})
    or a point pattern
    on a linear network (object of class \code{"lpp"}).
  }
}
\details{
  This function detects and repairs inconsistencies
  in the internal data of \code{X}. Currently it does the following:
  \itemize{
    \item checks that different ways of calculating the number of
    edges give the same answer
    \item removes any duplicated edges of the network
    \item ensures that each edge is recorded as a pair of vertex indices
    \code{(from, to)} with \code{from < to}.
  }
}
\value{
  An object of the same kind as \code{X}.
}
\author{
  \adrian.
}
\seealso{
  \code{\link{thinNetwork}}
}
\keyword{spatial}
\keyword{manip}
