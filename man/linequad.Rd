\name{linequad}
\alias{linequad}
\title{
  Quadrature Scheme on a Linear Network
}
\description{
  Generates a quadrature scheme (an object of class \code{"quad"})
  on a linear network.
}
\usage{
linequad(X, Y, \dots, eps = NULL, nd = 1000, random = FALSE)
}
\arguments{
  \item{X}{
    Data points. An object of class \code{"lpp"} or \code{"ppp"}.
  }
  \item{Y}{
    Line segments on which the points of \code{X} lie.
    An object of class \code{"psp"}.
    Required only when \code{X} is a \code{"ppp"} object.
  }
  \item{\dots}{
    Ignored.
  }
  \item{eps}{
    Optional. Spacing between successive dummy points along each
    segment. (This is the maximum spacing; some spacings will be shorter.)
  }
  \item{nd}{
    Optional. Total number of dummy points to be generated.
    (Actual number may be larger.)
  }
  \item{random}{
    Logical value indicating whether the sequence of dummy points
    should start at a randomly-chosen position along each segment.
  }
}
\details{
  This command generates a quadrature scheme (object of class
  \code{"quad"}) from a pattern of points on a linear network.

  Normally the user does not need to call \code{linequad} explicitly.
  It is invoked by \pkg{spatstat} functions when needed.
  A quadrature scheme is required by \code{\link{lppm}}
  in order to fit point process models to point pattern data on a linear
  network. A quadrature scheme is also used by \code{\link{rhohat.lpp}}
  and other functions.

  In order to create the quadrature scheme, dummy points are placed
  along each line segment of the network. The dummy points are 
  evenly-spaced with spacing \code{eps}. The default is
  \code{eps = totlen/nd} where \code{totlen} is the total length of
  all line segments in the network.

  Every line segment of the network will contain at least one dummy
  point. Consequently the actual number of dummy points generated
  will typically be greater than \code{nd},
  especially when \code{nd} is small.
  If \code{eps} is specified, the
  number of dummy points will be greater than \code{totlen/eps},
  especially when \code{eps} is large.
}
\value{
  A quadrature scheme (object of class \code{"quad"}).
}
\author{
  \adrian,
  Greg McSwiggan and Suman Rakshit.
}
\seealso{
  \code{\link{lppm}}
}
\keyword{datagen}
\keyword{spatial}
