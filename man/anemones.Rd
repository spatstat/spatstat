\name{anemones}
\alias{anemones}
\docType{data}
\title{
  Beadlet Anemones Data
}
\description{
  These data give the spatial locations and diameters
  of sea anemones (beadlet anemone
  \emph{Actinia equina}) in a sample plot 
  on the north face of a boulder,
  well above low tide level, at Quiberon (Bretagne, France) in May 1976.

  The data were originally described and discussed by Kooijman (1979a).
  Kooijman (1979b) shows a hand-drawn plot of the original data.
  The data are discussed  by Upton and Fingleton (1985)
  as Example 1.8 on pages 64--67. 

  The \code{anemones} dataset is taken directly from Table 1.11
  of Upton and Fingleton (1985). The coordinates and
  diameters are integer multiples of an idiosyncratic unit of length.
  The boundary is a rectangle 280 by 180 units.
}
\section{Units}{
  There is some confusion about the correct physical scale for these
  data. According to Upton and Fingleton (1985), one unit in the dataset 
  is approximately 0.475 cm. According to Kooijman (1979a, 1979b)
  and also quoted by Upton and Fingleton (1985), the
  physical size of the sample plot was 14.5 by 9.75 cm. 
  However if the data are plotted at this scale, they are too small for
  a rectangle of this size, and the appearance of the plot
  does not match the original hand-drawn plot in Kooijman (1979b).
  To avoid confusion, we have not assigned a unit scale to this
  dataset.
} 
\format{
  \code{anemones} is an object of class \code{"ppp"}
  representing the point pattern of anemone locations.
  It is a marked point pattern with numeric marks representing
  anemone diameter.
  See \code{\link{ppp.object}} for details of the format.
}
\usage{data(anemones)}
\examples{
  data(anemones)
  # plot diameters on same scale as x, y
  plot(anemones, markscale=0.5)
}
\source{
  Table 1.11 on pages 62--63 of Upton and Fingleton (1985),
  who acknowledge Kooijman (1979a) as the source.
}
\references{
 Kooijman, S.A.L.M. (1979a)
 The description of point patterns. 
 In \emph{Spatial and temporal analysis in ecology} (ed. R.M. Cormack
 and J.K. Ord), International Cooperative Publishing House, 
 Fairland, Maryland, USA. Pages 305--332.

 Kooijman, S.A.L.M. (1979b)
 Inference about dispersal patterns. 
 \emph{Acta Biotheoretica} \bold{28}, 149--189.

  Upton, G.J.G. and Fingleton, B. (1985)
  \emph{Spatial data analysis by example}. Volume 1: Point pattern
  and quantitative data. John Wiley and Sons, Chichester.
}
\keyword{datasets}
\keyword{spatial}
