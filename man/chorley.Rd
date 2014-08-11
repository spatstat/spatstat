\name{chorley}
\alias{chorley}
\alias{chorley.extra}
\docType{data}
\title{Chorley-Ribble Cancer Data}
\description{
  Spatial locations of cases of cancer of the larynx
  and cancer of the lung, and the location of a disused industrial
  incinerator. A marked point pattern.
} 
\format{
  The dataset \code{chorley} is
  an object of class \code{"ppp"}
  representing a marked point pattern.
  Entries include
  \tabular{ll}{
    \code{x} \tab Cartesian \eqn{x}-coordinate of home address \cr
    \code{y} \tab Cartesian \eqn{y}-coordinate of home address \cr
    \code{marks} \tab factor with levels \code{larynx} and \code{lung} \cr
    \tab indicating whether this is a case of cancer of the larynx\cr
    \tab or cancer of the lung.
  }
  See \code{\link{ppp.object}} for details of the format.

  The dataset \code{chorley.extra} is a list with two components.
  The first component \code{plotit} is a function which will 
  plot the data in a sensible fashion. The second
  component \code{incin} is a list with entries \code{x} and \code{y}
  giving the location of the industrial incinerator.

  Coordinates are given in kilometres,
  and the resolution is 100 metres (0.1 km)
}
\usage{data(chorley)}
\examples{
    chorley
    summary(chorley)
    chorley.extra$plotit()
}
\source{
  Coordinates of cases were provided by the
  Chorley and South Ribble Health Authority, and were
  kindly supplied by Professor Peter Diggle.
  Region boundary was digitised by Adrian Baddeley, 2005, from
  a photograph of an Ordnance Survey map.
}
\section{Notes}{
  The data give the precise domicile addresses of new cases 
  of cancer of the larynx (58 cases)
  and cancer of the lung (978 cases),
  recorded in the Chorley and South Ribble Health Authority
  of Lancashire (England) between 1974 and 1983.
  The supplementary data give the location of a disused industrial
  incinerator.

  The data were first presented and analysed by Diggle (1990).
  They have subsequently been analysed by Diggle and Rowlingson (1994)
  and Baddeley et al. (2005).
  
  The aim is to assess evidence for an increase in
  the incidence of cancer of the larynx in the vicinity of the
  now-disused industrial incinerator. The lung cancer cases serve as a
  surrogate for the spatially-varying density of the susceptible
  population.

  The data are represented as a marked point pattern,
  with the points giving the spatial location of each individual's home address
  and the marks identifying whether each point is a case of
  laryngeal cancer or lung cancer.

  Coordinates are in kilometres, and the resolution is
  100 metres (0.1 km).

  The dataset \code{chorley} has a polygonal window with 132 edges
  which closely approximates the boundary of the Chorley and South
  Ribble Health Authority. 

  Note that, due to the rounding of spatial coordinates,
  the data contain duplicated points (two points at the
  same location). To determine which points are duplicates,
  use \code{\link{duplicated.ppp}}.
  To remove the duplication, use \code{\link{unique.ppp}}.
}
\references{
  Baddeley, A., Turner, R., \ifelse{latex}{\out{M\o ller}}{Moller}, J. and Hazelton, M. (2005)
  Residual analysis for spatial point processes.
  \emph{Journal of the Royal Statistical Society, Series B}
  \bold{67}, 617--666.

  Diggle, P. (1990) A point process modelling approach to  
  raised incidence of a rare phenomenon in the vicinity
  of a prespecified point.
  \emph{Journal of the Royal Statistical Soc. Series A} 
  \bold{153}, 349-362.

  Diggle, P. and Rowlingson, B. (1994) A conditional approach
  to point process modelling of elevated risk.
  \emph{Journal of the Royal Statistical Soc. Series A} 
  \bold{157}, 433-440.
   
} 
\keyword{datasets}
\keyword{spatial}

