\name{bei}
\alias{bei}
\alias{bei.extra}
\docType{data}
\title{Tropical rain forest trees}
\description{
  A point pattern giving the locations of 3605
  trees in a tropical rain forest.
  Accompanied by covariate data giving the elevation (altitude)
  and slope of elevation in the study region.
} 
\format{
  \code{bei} is an object of class \code{"ppp"}
  representing the point pattern of tree locations.
  See \code{\link{ppp.object}} for details of the format.

  \code{bei.extra} is a list containing
  two pixel images, \code{elev} (elevation in metres) and
  \code{grad} (norm of elevation gradient). These pixel images are objects
  of class \code{"im"}, see \code{\link{im.object}}.
}
\usage{data(bei)}
\source{
  Hubbell and Foster (1983), Condit, Hubbell and Foster (1996)
  and Condit (1998).
  Data files kindly supplied by Rasmus Waagepetersen.
  The data were collected in the forest dynamics plot of Barro Colorado
  Island. The study was made possible through the generous support of
  the U.S. National Science Foundation, the John D. and Catherine
  T. MacArthur Foundation, and the Smithsonian Tropical Research Institute.
}
\section{Notes}{
  The dataset \code{bei} gives the positions of 3605 trees
  of the species \emph{Beilschmiedia pendula} (Lauraceae)
  in a 1000 by 500 metre rectangular sampling region
  in the tropical rainforest of Barro Colorado Island.

  The accompanying dataset  \code{bei.extra} gives information
  about the altitude (elevation) in the study region. It is a list
  containing two pixel images, \code{elev} (elevation in metres) and
  \code{grad} (norm of elevation gradient).

  These data are part of a much larger dataset containing the positions of
  hundreds of thousands of trees belong to thousands of species;
  see Hubbell and Foster (1983), Condit, Hubbell and Foster (1996)
  and Condit (1998).

  The present data were analysed by Moller and Waagepetersen (2007).
}
\references{
  Condit, R. (1998) \emph{Tropical Forest Census Plots}.
  Springer-Verlag, Berlin and R.G. Landes Company, Georgetown, Texas.

  Condit, R., Hubbell, S.P and Foster, R.B. (1996)
  Changes in tree species abundance in a neotropical forest: impact of
  climate change. \emph{Journal of Tropical Ecology} \bold{12},
  231--256.

  Hubbell, S.P and Foster, R.B. (1983)
  Diversity of canopy trees in a neotropical forest and implications for
  conservation. In: \emph{Tropical Rain Forest: Ecology and Management}
  (eds. S.L. Sutton, T.C. Whitmore and A.C. Chadwick),
  Blackwell Scientific Publications, Oxford, 25--41.

  Moller, J. and Waagepetersen, R.P. (2007)
  Modern spatial point process modelling and inference (with discussion).
  \emph{Scandinavian Journal of Statistics} \bold{34}, 643--711.
} 
\keyword{datasets}
\keyword{spatial}

