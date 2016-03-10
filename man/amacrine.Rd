\name{amacrine}
\alias{amacrine}
\docType{data}
\title{Hughes' Amacrine Cell Data}
\description{
Austin Hughes' data: a point pattern 
of displaced amacrine cells in the retina of a rabbit.
A marked point pattern.
} 
\format{
  An object of class \code{"ppp"}
  representing the point pattern of cell locations.
  Entries include
  \tabular{ll}{
    \code{x} \tab Cartesian \eqn{x}-coordinate of cell \cr
    \code{y} \tab Cartesian \eqn{y}-coordinate of cell \cr
    \code{marks} \tab factor with levels \code{off} and \code{on} \cr
                 \tab indicating ``off'' and ``on'' cells
    }
  See \code{\link{ppp.object}} for details of the format.
}
\usage{data(amacrine)}
\source{Peter Diggle, personal communication}
\section{Notes}{
Austin Hughes' data: a point pattern 
of displaced amacrine cells in the retina of a rabbit.
152 ``on'' cells and 142 ``off'' cells in a rectangular sampling frame.

The true dimensions of the rectangle are 1060 by 662 microns.
The coordinates here are scaled to a rectangle of height 1 and width
\eqn{1060/662 = 1.601} so the unit of measurement is approximately 662 microns.

The data were analysed by Diggle (1986).
}
\references{
Diggle, P. J. (1986).
Displaced amacrine cells in the retina of a
rabbit: analysis of a bivariate spatial point pattern. 
\emph{J. Neurosci. Meth.} \bold{18}, 115--125.
} 
\keyword{datasets}
\keyword{spatial}

