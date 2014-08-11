\name{hamster}
\alias{hamster}
\docType{data}
\title{Aherne's hamster tumour data} 
\description{
  Point pattern of cell nuclei in hamster kidney, each nucleus classified as
  either `dividing' or `pyknotic'.
  A multitype point pattern.
} 
\format{
  An object of class \code{"ppp"}
  representing the point pattern of cell locations.
  Entries include
  \tabular{ll}{
    \code{x} \tab Cartesian \eqn{x}-coordinate of cell \cr
    \code{y} \tab Cartesian \eqn{y}-coordinate of cell \cr
    \code{marks} \tab factor with levels \code{"dividing"}
    and \code{"pyknotic"}.
  }
  See \code{\link{ppp.object}} for details of the format.
}
\usage{data(hamster)}
\source{Dr W. A. Aherne, Department of Pathology, University of
  Newcastle-upon-Tyne, UK. Data supplied by Prof. Peter Diggle}
\section{Notes}{
  These data were presented and analysed by Diggle (1983, section 7.3).

  The data give the positions of the centres of the nuclei of certain
  cells in a histological
  section of tissue from a laboratory-induced metastasising lymphoma
  in the kidney of a hamster.

  The nuclei are classified as either
  "pyknotic" (corresponding to dying cells) or "dividing" (corresponding
  to cells arrested in metaphase, i.e. in the act of dividing). The
  background void is occupied by unrecorded, interphase cells in
  relatively large numbers. 

  The sampling window is a square, originally about 0.25 mm square
  in real units, which has been rescaled to the unit square.
}
\references{
  Diggle, P.J. (1983)
  \emph{Statistical analysis of spatial point patterns}.
  Academic Press.
}
\keyword{datasets}
\keyword{spatial}
