\name{sporophores}
\alias{sporophores}
\docType{data}
\title{
  Sporophores Data
}
\description{
  Spatial pattern of sporophores of three species of fungi around a tree.
}
\usage{data(sporophores)}
\format{
  A multitype spatial point pattern (an object of class \code{"ppp"}
  with factor-valued marks indicating the species).
  Spatial coordinates are given in centimetres.
  Levels of the species variable are
  \code{"L laccata"}, \code{"L pubescens"} and \code{"Hebloma spp"}.
}
\details{
  Ford, Mason and Pelham (1980) studied the spatial locations
  of sporophores of three species of mycorrhizal fungi
  distributed around a young birch tree in agricultural soil.
  The dataset given here is the spatial pattern
  in the fifth year after the tree was planted.
  The species are 
  \emph{Laccaria laccata}, \emph{Lactarius pubescens}
  and \emph{Hebloma} spp.
}
\source{
  Data generously provided by Dr E.D. Ford.
  Please cite Ford et al (1980) in any use of these data.
}
\references{
  Ford, E.D., Mason, P.A. and Pelham, J. (1980) 
  Spatial patterns of sporophore distribution around a
  young birch tree in three successive years.
  \emph{Transactions of the British Mycological Society}
  \bold{75}, 287--296.
}
\examples{
## reproduce Fig 1 in Ford et al (1980)
plot(sporophores, chars=c(16,1,2), cex=0.6, leg.args=list(cex=1.1))
points(0,0,pch=16, cex=2)
text(15,8,"Tree", cex=0.75)
}
\keyword{datasets}
