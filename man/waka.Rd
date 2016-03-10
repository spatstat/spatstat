\name{waka}
\alias{waka}
\docType{data}
\title{
  Trees in Waka national park
}
\description{
  This dataset is a spatial point pattern of trees 
  recorded at Waka National Park, Gabon.
  See Balinga et al (2006).

  The dataset \code{waka} is a point pattern
  (object of class \code{"ppp"}) containing the spatial coordinates
  of each tree, marked by the tree diameter at breast height
  \code{dbh}.
  The survey region is a 100 by 100 metre square.
  Coordinates are given in metres, while the \code{dbh} is in centimetres.
}
\usage{data(waka)}
\examples{
data(waka)
plot(waka, markscale=0.01)
title(sub="Tree diameters to scale")
plot(waka, markscale=0.04)
title(sub="Tree diameters 4x scale")
}
\source{
  Nicolas Picard
}
\references{
  Balinga, M., Sunderland, T., Walters, G., Issemb{\'e}, Y., Asaha, S. 
  and Fombod, E. (2006)
  \emph{A vegetation assessment of the Waka national park, Gabon.}
  Herbier National du Gabon, LBG, MBG, WCS, FRP and
  Simthsonian Institution, Libreville, Gabon. CARPE Report, 154 pp.
  \url{http://carpe.umd.edu/}

  Picard, N., Bar-Hen, A., Mortier, F. and Chadoeuf, J. (2009)
  The multi-scale marked area-interaction point process: a model for
  the spatial pattern of trees.
  \emph{Scandinavian Journal of Statistics} \bold{36} 23--41
}
\keyword{datasets}
\keyword{spatial}
