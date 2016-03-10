\name{ponderosa}
\alias{ponderosa}
\alias{ponderosa.extra}
\docType{data}
\title{
  Ponderosa Pine Tree Point Pattern
}
\description{
  The data record the locations of 108 Ponderosa Pine
  (\emph{Pinus ponderosa}) trees in a 120 metre square region in the Klamath
  National Forest in northern California, published as Figure 2 of
  Getis and Franklin (1987).

  Franklin et al. (1985) determined the locations of approximately 5000
  trees from United States Forest Service aerial photographs and
  digitised them for analysis. Getis and Franklin (1987) selected a 120
  metre square subregion that appeared to exhibit clustering. This subregion
  is the \code{ponderosa} dataset.

  In principle these data are equivalent to
  Figure 2 of Getis and Franklin (1987) but they are not exactly
  identical; some of the spatial locations appear to be slightly
  perturbed.

  The data points identified as A, B, C on Figure 2 of Getis and Franklin
  (1987) correspond to points numbered 42, 7 and 77 in the dataset
  respectively.
} 
\format{
  Typing \code{data(ponderosa)} gives access to two objects,
  \code{ponderosa} and \code{ponderosa.extra}.
  
  The dataset \code{ponderosa} is a spatial point pattern 
  (object of class \code{"ppp"})
  representing the point pattern of tree positions.
  See \code{\link{ppp.object}} for details of the format.

  The dataset \code{ponderosa.extra} is a list containing supplementary
  data. The entry \code{id} contains the index numbers of the
  three special points A, B, C in the point pattern. The entry
  \code{plotit} is a function that can be called to produce a nice plot
  of the point pattern.
}
\usage{data(ponderosa)}
\source{Prof. Janet Franklin, University of California, Santa Barbara}
\examples{
   data(ponderosa)
   ponderosa.extra$plotit()
}
\references{
  Franklin, J., Michaelsen, J. and Strahler, A.H. (1985)
  Spatial analysis of density dependent pattern in coniferous forest
  stands.
  \emph{Vegetatio} \bold{64}, 29--36.
  
  Getis, A. and Franklin, J. (1987)
  Second-order neighbourhood analysis of mapped point patterns.
  \emph{Ecology} \bold{68}, 473--477.
}
\keyword{datasets}
\keyword{spatial}
