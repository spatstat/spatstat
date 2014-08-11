\name{mucosa}
\alias{mucosa}
\alias{mucosa.subwin}
\docType{data}
\title{
  Cells in Gastric Mucosa
}
\description{
  A bivariate inhomogeneous point pattern, giving the locations of
  the centres of two types of cells in a cross-section of the
  gastric mucosa of a rat.
}
\usage{data(mucosa)}
\format{
  An object of class \code{"ppp"}, see \code{\link{ppp.object}}.
  This is a multitype point pattern with two types of points,
  \code{ECL} and \code{other}.
}
\details{
  This point pattern dataset gives the locations of 
  cell centres in a cross-section of the gastric mucosa (mucous membrane of the
  stomach) of a rat. The rectangular observation window has been scaled
  to unit width. The lower edge of the window is closest to the outside
  of the stomach.

  The cells are classified into two types: \emph{ECL cells}
  (enterochromaffin-like cells) and \emph{other} cells. There are 86 ECL cells
  and 807 other cells in the dataset. 
  ECL cells are a type of neuroendocrine cell which 
  synthesize and secrete histamine.
  One hypothesis of interest is whether the 
  spatially-varying intensities of ECL cells and other cells are
  proportional.

  The data were originally collected by Dr Thomas Berntsen.
  The data were discussed and analysed in \ifelse{latex}{\out{M\o ller}}{Moller} and Waagepetersen (2004,
  pp. 2, 169).

  The associated object \code{mucosa.subwin} is the smaller window
  to which the data were restricted for analysis by \ifelse{latex}{\out{M\o ller}}{Moller} and Waagepetersen.
}
\source{
  Dr Thomas Berntsen and Prof Rasmus Waagepetersen.
}
\references{
  \ifelse{latex}{\out{M\o ller}}{Moller}, J. and Waagepetersen, R. (2004).
  \emph{Statistical Inference and Simulation for Spatial Point Processes}.
  Chapman and Hall/CRC.
}
\examples{
  plot(mucosa, chars=c(1,3), cols=c("red", "green"))
  plot(mucosa.subwin, add=TRUE, lty=3)
}
\keyword{datasets}
