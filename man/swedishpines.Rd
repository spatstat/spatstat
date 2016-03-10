\name{swedishpines}
\alias{swedishpines}
\docType{data}
\title{
  Swedish Pines Point Pattern
}
\description{
  The data give the locations of pine saplings
  in a Swedish forest. 
} 
\format{
  An object of class \code{"ppp"}
  representing the point pattern of tree locations
  in a rectangular plot 9.6 by 10 metres.

  Cartesian coordinates are given in decimetres (multiples of 0.1 metre)
  rounded to the nearest decimetre.
  Type \code{\link{rescale}(swedishpines)} to get an equivalent dataset
  where the coordinates are expressed in metres.
  
  See \code{\link{ppp.object}} for details of the format of a
  point pattern object.
}
\usage{data(swedishpines)}
\source{Strand (1972), Ripley (1981)}
\section{Note}{
  For previous analyses see Ripley (1981, pp. 172-175),
  Venables and Ripley (1997, p. 483),
  Baddeley and Turner (2000).
}
\references{
  Baddeley, A. and Turner, R. (2000)
  Practical maximum pseudolikelihood for spatial point patterns.
  \emph{Australian and New Zealand Journal of Statistics}
  \bold{42}, 283--322.
 
  Ripley, B.D. (1981) 
  \emph{Spatial statistics}.
  John Wiley and Sons.

  Strand, L. (1972). 
  A model for stand growth.
  \emph{IUFRO Third Conference Advisory Group of Forest Statisticians},
  INRA, Institut National de la Recherche Agronomique, Paris.
  Pages 207--216.

  Venables, W.N. and Ripley, B.D. (1997)
  \emph{Modern applied statistics with S-PLUS}.
  Second edition. Springer Verlag. 
}
\keyword{datasets}
\keyword{spatial}
