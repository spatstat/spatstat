\name{hyytiala}
\alias{hyytiala}
\docType{data}
\title{
  Scots pines and other trees at Hyytiala
}
\description{
  This dataset is a spatial point pattern of trees recorded at
  \ifelse{latex}{\out{Hyyti\"{a}l\"{a}}}{Hyytiala}, Finland. 
  The majority of the trees are Scots pines.
  See Kokkila et al (2002).

  The dataset \code{hyytiala} is a point pattern
  (object of class \code{"ppp"}) containing the spatial coordinates
  of each tree, marked by species (a factor with levels \code{aspen},
  \code{birch}, \code{pine} and \code{rowan}).
  The survey region is a 20 by 20 metre square.
  Coordinates are given in metres.
}
\usage{data(hyytiala)}
\examples{
data(hyytiala)
plot(hyytiala, cols=2:5)
}
\source{
  Nicolas Picard
}
\references{
  Kokkila, T., \ifelse{latex}{\out{M{\"a}kel{\"a}}}{Makela}, A. and
  Nikinmaa E. (2002)
  A method for generating stand structures using Gibbs marked point
  process. 
  \emph{Silva Fennica} \bold{36} 265--277.

  Picard, N, Bar-Hen, A., Mortier, F. and Chadoeuf, J. (2009)
  The multi-scale marked area-interaction point process: a model for
  the spatial pattern of trees.
  \emph{Scandinavian Journal of Statistics} \bold{36} 23--41
}
\keyword{datasets}
\keyword{spatial}
