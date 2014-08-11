\encoding{latin1}
\name{urkiola}
\alias{urkiola}
\docType{data}
\title{Urkiola Woods Point Pattern}
\description{
  Locations of birch (\emph{Betula celtiberica}) and
  oak (\emph{Quercus robur}) trees  in a secondary wood in
  Urkiola Natural Park (Basque country, northern Spain). 
  They are part of a more extensive dataset collected
  and analysed by Laskurain (2008). The coordinates of the trees
  are given in meters. 
}
\usage{data(urkiola)}
\format{
  An object of class \code{"ppp"} representing the point pattern of
  tree locations. Entries include
  \describe{
    \item{x}{Cartesian x-coordinate of tree}
    \item{y}{Cartesian y-coordinate of tree }
    \item{marks}{factor indicating species of each tree}
  }
  The levels of \code{marks} are \code{birch} and \code{oak}.
  See \code{\link{ppp.object}}  for details of the format of a ppp object. 
}
\source{N.A. Laskurain. Kindly formatted and communicated by M. de la Cruz Rot}
\references{
  Laskurain, N. A. (2008)
  \emph{\enc{Dinámica}{Dinamica} espacio-temporal de un bosque
    secundario en el Parque Natural de Urkiola (Bizkaia).} 
  Tesis Doctoral. Universidad del \enc{País}{Pais} Vasco /Euskal Herriko
  Unibertsitatea.
}
\keyword{datasets}
