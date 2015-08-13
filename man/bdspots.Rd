\name{bdspots}
\alias{bdspots}
\docType{data}
\title{
  Breakdown Spots in Microelectronic Materials
}
\description{
  A list of three point patterns, each giving the locations of
  electrical breakdown spots on a circular electrode in
  a microelectronic capacitor.
}
\usage{data(bdspots)}
\format{
  A list (of class \code{"listof"}) of three spatial point patterns,
  each representing the spatial locations of breakdown spots on an
  electrode. The three electrodes are circular discs, of radii
  169, 282 and 423 microns respectively. Spatial coordinates are
  given in microns.  
}
\details{
  The application of successive voltage sweeps to the metal gate electrode
  of a microelectronic capacitor generates multiple breakdown spots
  on the electrode.
  The spatial distribution of these breakdown spots
  in MIM (metal-insulator-metal) and MIS (metal-insulator-semiconductor)
  structures was observed and analysed by 
  Miranda et al (2010, 2013) and Saura et al (2013a, 2013b, 2014).

  The data given here are the breakdown spot patterns for three circular
  electrodes of different radii, 169, 282 and 423 microns respectively,
  in MIM structures analysed in Saura et al (2013a).
}
\source{
  Professor Enrique Miranda, Departament d'Enginyeria
  \ifelse{latex}{\out{Electr{\`o}nica}}{Electronica},
  Escola d'Enginyeria, Universitat
  \ifelse{latex}{\out{Aut{\`o}noma}}{Autonoma}
  de Barcelona, Barcelona, Spain.
}
\references{
  Miranda, E. and O'Connor, E. and Hurley, P.K. (2010)
  Simulation of the breakdown spots spatial
  distribution in high-\emph{K} dielectrics and model
  validation using the \pkg{spatstat} package for \emph{R} language.
  \emph{ECS Transactions} \bold{33} (3) 557--562.

  Miranda, E.,
  \ifelse{latex}{\out{Jim{\'e}nez}}{Jimenez}, D.,
  \ifelse{latex}{\out{Su{\~n}{\'e}}}{Sune}, J.,
  O'Connor, E.,
  Monaghan, S.,
  Povey, I.,
  Cherkaoui, K. and Hurley, P. K. (2013)
  Nonhomogeneous spatial distribution of filamentary leakage current
  paths in circular area Pt/HfO2/Pt capacitors.
  \emph{J. Vac. Sci. Technol. B} \bold{31}, 01A107.

  Saura, X.,
  \ifelse{latex}{\out{Su{\~n}{\'e}}}{Sune}, J.,
  Monaghan, S., Hurley, P.K. and Miranda, E. (2013a)
  Analysis of the breakdown spot spatial distribution in Pt/HfO2/Pt
  capacitors using nearest neighbor statistics.
  \emph{J. Appl. Phys.} \bold{114}, 154112.

  Saura, X., Moix, D.,
  \ifelse{latex}{\out{Su{\~n}{\'e}}}{Sune}, J.,
  Hurley, P.K. and Miranda, E. (2013b)
  Direct observation of the generation of breakdown spots in MIM
  structures under constant voltage stress.
  \emph{Microelectronics Reliability} \bold{53}, 1257--1260.

  Saura, X.,
  Monaghan, S.,
  Hurley, P.K.,
  \ifelse{latex}{\out{Su{\~n}{\'e}}}{Sune}, J.
  and Miranda, E. (2014)
  Failure analysis of MIM and MIS structures using point-to-event
  distance and angular probability distributions. 
  \emph{IEEE Transactions on Devices and Materials Reliability} 
  \bold{14} (4) 1080--1090.
}
\examples{
data(bdspots)
plot(bdspots, equal.scales=TRUE)
}
\keyword{datasets}
\keyword{spatial}
