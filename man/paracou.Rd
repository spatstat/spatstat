\name{paracou}
\alias{paracou}
\docType{data}
\title{
  Kimboto trees at Paracou, French Guiana
}
\description{
  This dataset is a point pattern of adult and juvenile Kimboto trees
  (\emph{Pradosia cochlearia} or \emph{P. ptychandra})
  recorded at Paracou in French Guiana.
  See Flores (2005).

  The dataset \code{paracou} is a point pattern
  (object of class \code{"ppp"}) containing the spatial coordinates
  of each tree, marked by age (a factor with levels \code{adult} and
  \code{juvenile}. The survey region is a rectangle
  approximately 400 by 525 metres. Coordinates are given in metres.

  Note that the data contain duplicated points (two points at the
  same location). To determine which points are duplicates,
  use \code{\link{duplicated.ppp}}.
  To remove the duplication, use \code{\link{unique.ppp}}.
}
\usage{data(paracou)}
\examples{
plot(paracou, cols=2:3, chars=c(16,3))
}
\source{
  Data kindly contributed by Olivier Flores.
  All data belong to CIRAD \url{http://www.cirad.fr}
  and UMR EcoFoG \url{http://www.ecofog.gf} and
  are included in \pkg{spatstat} with permission.
  Original data sources:
  juvenile and some adult trees collected by Flores (2005);
  adult tree data sourced from CIRAD Paracou experimental plots dataset
  (2003 campaign).
}
\references{
  Flores, O. (2005)
  \emph{\ifelse{latex}{\out{D{\'e}terminisme de la
  r{\'e}g{\'e}n{\'e}ration chez quinze esp{\`e}ces d'arbres tropicaux
  en for{\^e}t guyanaise: les effets de l'environnement et de la
  limitation par la dispersion.}}{ Determinisme de la regeneration
  chez quinze especes d'arbres tropicaux en foret guyanaise: les
  effets de l'environnement et de la limitation par la dispersion.}}
  PhD Thesis, University of Montpellier 2, Montpellier, France.

  Picard, N, Bar-Hen, A., Mortier, F. and Chadoeuf, J. (2009)
  The multi-scale marked area-interaction point process: a model for
  the spatial pattern of trees.
  \emph{Scandinavian Journal of Statistics} \bold{36} 23--41
}
\keyword{datasets}
\keyword{spatial}
