\name{chicago}
\alias{chicago}
\docType{data}
\title{
  Chicago Street Crime Data
}
\description{
  This dataset is a record of street crimes
  reported in the period 25 April to 8 May 2002,
  in an area of Chicago (Illinois, USA)
  close to the University of Chicago.
  The original street crime map was published in the
  Chicago Weekly News in 2002. 
  
  The data give the spatial location of each crime,
  and the type of crime. The type labels are interpreted as follows:
  \tabular{ll}{
    \code{assault} & battery/assault \cr
    \code{burglary} & burglary \cr
    \code{cartheft} & motor vehicle theft \cr
    \code{damage} & criminal damage \cr
    \code{robbery} & robbery \cr
    \code{theft} & theft \cr
    \code{trespass} & criminal trespass
  }
  All crimes occurred on or near a street. The data give the
  coordinates of all streets in the survey area, and their connectivity.

  The dataset \code{chicago} is an object of class \code{"lpp"}
  representing a point pattern on a linear network.
  See \code{\link{lpp}} for further information on the format.

  These data were published and analysed in
  Ang, Baddeley and Nair (2012).
} 
\format{
  Object of class \code{"lpp"}. 
  See \code{\link{lpp}}.
}
\usage{data(chicago)}
\examples{
data(chicago)
plot(chicago)
plot(as.linnet(chicago), main="Chicago Street Crimes",col="green")
plot(as.ppp(chicago), add=TRUE, col="red", chars=c(16,2,22,17,24,15,6))
}
\source{
  Chicago Weekly News, 2002.
  Manually digitised by Adrian Baddeley.
}
\references{
  Ang, Q.W. (2010)
  \emph{Statistical methodology for events on a network}.
  Master's thesis, School of Mathematics and Statistics, University of
  Western Australia.
  
  Ang, Q.W., Baddeley, A. and Nair, G. (2012)
  Geometrically corrected second-order analysis of 
  events on a linear network, with applications to
  ecology and criminology.
  \emph{Scandinavian Journal of Statistics} \bold{39}, 591--617.

  Chicago Weekly News website: \url{www.chicagoweeklynews.com}
}
\keyword{datasets}
\keyword{spatial}
