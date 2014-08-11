\name{longleaf}
\alias{longleaf}
\docType{data}
\title{
  Longleaf Pines Point Pattern
}
\description{
  Locations and sizes of Longleaf pine trees.
  A marked point pattern.

  The data record the locations and diameters of 
  584 Longleaf pine (\emph{Pinus palustris}) trees 
  in a 200 x 200 metre region in southern Georgia (USA).
  They were collected and analysed by Platt, Evans and Rathbun (1988).
 
  This is a marked point pattern; the mark associated with a tree is its
  diameter at breast height (\code{dbh}), a convenient measure of its size. 
  Several analyses have considered only the ``adult'' trees which
  are conventionally defined as those trees with \code{dbh}
  greater than or equal to 30 cm.
 
  The pattern is regarded as spatially inhomogeneous.
} 
\format{
  An object of class \code{"ppp"}
  representing the point pattern of tree locations.
  Entries include
  \tabular{ll}{
    \code{x} \tab Cartesian \eqn{x}-coordinate of tree \cr
    \code{y} \tab Cartesian \eqn{y}-coordinate of tree \cr
    \code{marks} \tab diameter at breast height, in centimetres.
  }
  See \code{\link{ppp.object}} for details of the format of a
  point pattern object.
}
\usage{data(longleaf)}
\examples{
    data(longleaf)
    plot(longleaf)
    plot(cut(longleaf, breaks=c(0,30,Inf), labels=c("Sapling","Adult")))
}
\source{Platt, Evans and Rathbun (1988)}
\references{
Platt, W. J., Evans, G. W. and Rathbun, S. L. (1988)
The population dynamics of a long-lived Conifer (Pinus palustris).
\emph{The American Naturalist} \bold{131}, 491--525.

Rathbun, S. L. and Cressie, N. (1994)
A space-time survival point process for a longleaf
pine forest in southern Georgia.
\emph{Journal of the American Statistical Association} \bold{89}, 1164--1173.
}
\keyword{datasets}
\keyword{spatial}
