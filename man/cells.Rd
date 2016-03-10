\name{cells}
\alias{cells}
\docType{data}
\title{
  Biological Cells Point Pattern
}
\description{
The data record the locations of the centres of 42 biological cells
observed under optical microscopy in a histological section. 
The microscope field-of-view has been rescaled to the unit square.

The data were recorded by F.H.C. Crick and B.D. Ripley,
and analysed in Ripley (1977, 1981) and Diggle (1983).
They are often used as a canonical example of an `ordered'
point pattern.
} 
\format{
  An object of class \code{"ppp"}
  representing the point pattern of cell centres.
  See \code{\link{ppp.object}} for details of the format.
}
\usage{data(cells)}
\source{Crick and Ripley, see Ripley (1977)}
\references{
  Diggle, P.J. (1983)
  \emph{Statistical analysis of spatial point patterns}.
  Academic Press.

  Ripley, B.D. (1977)
  Modelling spatial patterns (with discussion).
  \emph{Journal of the Royal Statistical Society, Series B}
  \bold{39}, 172--212.

  Ripley, B.D. (1981) 
  \emph{Spatial statistics}.
  John Wiley and Sons.
}
\keyword{datasets}
\keyword{spatial}
