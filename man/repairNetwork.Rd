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
  in the internal data of \code{X}. Currently it detects and removes
  \itemize{
    \item duplication of edges of the network (the presence of two
    edges between the same pair of vertices).
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
