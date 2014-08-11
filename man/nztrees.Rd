\name{nztrees}
\alias{nztrees}
\docType{data}
\title{
   New Zealand Trees Point Pattern
}
\description{
  The data give the locations of trees in a forest plot.

  They were collected by Mark and Esler (1970) and
  were extracted and analysed by Ripley (1981, pp. 169-175). 
  They represent the 
  positions of 86 trees in a forest plot approximately 140 by 85 feet.
  
  Ripley discarded from his analysis the eight trees at the right-hand edge of
  the plot (which appear to be part of a planted border)
  and trimmed the window by a 5-foot margin accordingly.
} 
\format{
  An object of class \code{"ppp"}
  representing the point pattern of tree locations.
  The Cartesian coordinates are in feet.

  See \code{\link{ppp.object}} for details of the format of a
  point pattern object.
}
\usage{data(nztrees)}
\source{Mark and Esler (1970), Ripley (1981).}
\section{Note}{
  To trim a 5-foot margin off the window, type
  \code{nzsub <- nztrees[ , owin(c(0,148),c(0,95)) ]}
}
\references{
  Ripley, B.D. (1981) 
  \emph{Spatial statistics}.
  John Wiley and Sons.

  Mark, A.F. and Esler, A.E. (1970)
  An assessment of the point-centred quarter method of plotless
  sampling in some New Zealand forests.
  \emph{Proceedings of the New Zealand Ecological Society}
  \bold{17}, 106--110.
}
\keyword{datasets}
\keyword{spatial}


 
