\name{murchison}
\alias{murchison}
\docType{data}
\title{
  Murchison gold deposits
}
\description{
  Data recording the spatial locations of gold deposits and associated
  geological features in the Murchison area of Western Australia.
  Extracted from a large scale (1:500,000) study of the
  Murchison area by the Geological Survey of Western Australia
  (Watkins and Hickman, 1990). The features recorded are
  \itemize{
    \item the locations of gold deposits;
    \item the locations of geological faults;
    \item the region that contains greenstone bedrock.
  }
  The study region is contained in a \eqn{330\times 400}{330 * 400} kilometre
  rectangle. At this scale, gold deposits are points, i.e. their spatial
  extent is negligible. 
  Gold deposits in this region occur only in greenstone bedrock.
  Geological faults can be observed reliably only within the same
  region. However, some faults have been extrapolated
  (by geological ``interpretation'') outside the greenstone boundary 
  from information observed in the greenstone region.

  These data were analysed by Foxall and Baddeley (2002) and Brown et al
  (2002); see also Groves et al (2000), Knox-Robinson and Groves (1997).
  The main aim is to predict the intensity of the
  point pattern of gold deposits from the more easily observable fault
  pattern.
}
\format{
  \code{murchison} is a list with the following entries:
  \describe{
    \item{gold}{a point pattern (object of class \code{"ppp"})
      representing the point pattern of gold deposits.
      See \code{\link{ppp.object}} for details of the format.
    }
    \item{faults}{a line segment pattern (object of class \code{"psp"})
      representing the geological faults.
      See \code{\link{psp.object}} for details of the format.
    }
    \item{greenstone}{the greenstone bedrock region.
      An object of class \code{"owin"}. Consists of multiple
      irregular polygons with holes.
    }
  }
  All coordinates are given in \bold{metres}.
}
\usage{
 data(murchison)
}
\examples{
  if(interactive()) {
  data(murchison)
  plot(murchison$greenstone, main="Murchison data", col="lightgreen")
  plot(murchison$gold, add=TRUE, pch="+",col="blue")
  plot(murchison$faults, add=TRUE, col="red")
  }
}
\source{
  Data were kindly provided by Dr Carl Knox-Robinson of the
  Department of Geology and Geophysics, University of Western Australia.
  Permission to use the data is granted by Dr Tim Griffin,
  Geological Survey of Western Australia and by Dr Knox-Robinson.
  \emph{Please make appropriate acknowledgement} to
  Watkins and Hickman (1990) and the Geological Survey of Western Australia.
}
\references{
  Brown, W.M., Gedeon, T.D., Baddeley, A.J. and Groves, D.I. (2002)
  Bivariate J-function and other graphical
  statistical methods help select the best predictor
  variables as inputs for a neural network method of
  mineral prospectivity mapping.
  In U. Bayer, H. Burger and W. Skala (eds.)
  \emph{IAMG 2002: 8th Annual Conference of the
    International Association for Mathematical Geology},
  Volume 1, 2002. International Association of Mathematical Geology.
  Pages 257--268.

  Foxall, R. and Baddeley, A. (2002)
  Nonparametric measures of association between a
  spatial point process and a random set, with
  geological applications. \emph{Applied Statistics} \bold{51}, 165--182.

  Groves, D.I., Goldfarb, R.J., Knox-Robinson, C.M., Ojala, J., Gardoll,
  S, Yun, G.Y. and Holyland, P. (2000)
  Late-kinematic timing of orogenic gold deposits and significance for
  computer-based exploration techniques with emphasis on the Yilgarn Block,
  Western Australia.
  \emph{Ore Geology Reviews}, \bold{17}, 1--38.

  Knox-Robinson, C.M. and Groves, D.I. (1997)
  Gold prospectivity mapping using a geographic information system
  (GIS), with examples from the Yilgarn Block of Western Australia.
  \emph{Chronique de la Recherche Miniere} \bold{529}, 127--138.

  Watkins, K.P. and Hickman, A.H. (1990)
  \emph{Geological evolution and mineralization of the Murchison Province,
    Western Australia}.
  Bulletin 137, Geological Survey of Western Australia. 267 pages.
  Published by Department of Mines, Western Australia, 1990.
  Available online from Department of Industry and Resources,
  State Government of Western Australia, \code{www.doir.wa.gov.au}
}
\keyword{datasets}
\keyword{spatial}
