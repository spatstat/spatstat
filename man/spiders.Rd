\name{spiders}
\alias{spiders}
\docType{data}
\title{
  Spider Webs on Mortar Lines of a Brick Wall
}
\description{
  Data recording the locations of small spider webs on
  the network of mortar lines of a brick wall.
}
\usage{data("spiders")}
\format{
  Object of class \code{"lpp"} representing a pattern of points
  on a linear network.
}
\details{
The data give the positions of 48 webs 
of the urban wall spider \emph{Oecobius navus} 
on the mortar lines of a brick wall, 
recorded by Voss (1999) and manually digitised by Mark Handcock. 
The mortar spaces provide the only opportunity for constructing webs
(Voss 1999; Voss et al 2007) so this is a pattern of points on a
network of lines.

The habitat preferences of this species were studied in detail
by Voss et al (2007). Questions of interest include evidence for
non-uniform density of webs and for interaction between nearby individuals.

Observations were made inside a square quadrat of side length 1.125 metres.
The original hand-drawn map was digitised manually by Mark
S. Handcock, and reformatted as a \code{spatstat} object by
Ang Qi Wei.

The dataset \code{spiders} is an object of class \code{"lpp"}
(point pattern on a linear network). Coordinates are given in millimetres.
The linear network has 156 vertices and a total length of 20.22 metres. 

\emph{Please cite Voss et al (2007) with any use of these data.} 
}
\source{
  Dr Sasha Voss. Coordinates manually recorded by M.S. Handcock
  and formatted by Q.W. Ang.
  
  \emph{Please cite Voss et al (2007) with any use of these data.} 
}
\references{
  Ang, Q.W. (2010)
  \emph{Statistical methodology for events on a network}.
  Master's thesis, School of Mathematics and Statistics, University of
  Western Australia.

  Voss, S. (1999)
  Habitat preferences and spatial dynamics of the
  urban wall spider: \emph{Oecobius annulipes} Lucas.
  Honours Thesis, Department of Zoology, University of Western
   Australia.	  

  Voss, S., Main, B.Y. and Dadour, I.R. (2007)
  Habitat preferences of the urban wall spider 
  \emph{Oecobius navus} (Araneae, Oecobiidae).
  \emph{Australian Journal of Entomology} \bold{46}, 261--268.
}
\examples{
plot(spiders, show.window=FALSE, pch=16)
}
\keyword{datasets}
