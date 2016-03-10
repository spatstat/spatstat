\name{lansing}
\alias{lansing}
\docType{data}
\title{
  Lansing Woods Point Pattern
}
\description{
  Locations and botanical classification of trees in Lansing Woods.

  The data come from an investigation of a 924 ft x 924 ft (19.6 acre)
  plot in Lansing Woods, Clinton County, Michigan USA
  by D.J. Gerrard. The data give the locations of 2251 trees and 
  their botanical classification (into hickories, maples, red oaks,
  white oaks, black oaks and miscellaneous trees).
  The original plot size (924 x 924 feet)
  has been rescaled to the unit square.
  
  Note that the data contain duplicated points (two points at the
  same location). To determine which points are duplicates,
  use \code{\link{duplicated.ppp}}.
  To remove the duplication, use \code{\link{unique.ppp}}.
} 
\format{
  An object of class \code{"ppp"}
  representing the point pattern of tree locations.
  Entries include
  \tabular{ll}{
    \code{x} \tab Cartesian \eqn{x}-coordinate of tree \cr
    \code{y} \tab Cartesian \eqn{y}-coordinate of tree \cr
    \code{marks} \tab factor with levels indicating species of
      each tree
  }
  The levels of \code{marks} are
  \code{blackoak},
  \code{hickory},
  \code{maple},
  \code{misc},
  \code{redoak} and
  \code{whiteoak}.
  See \code{\link{ppp.object}} for details of the format of a
  point pattern object.
}
\usage{data(lansing)}
\examples{
     data(lansing)
     plot(lansing)
     summary(lansing)
     plot(split(lansing))
     plot(split(lansing)$maple)
}
\references{
Besag, J. (1978)
Some methods of statistical analysis for spatial data.
\emph{Bull. Internat. Statist. Inst.} \bold{44}, 77--92.

Cox, T.F. (1976)
The robust estimation of the density of a forest
stand using a new conditioned distance method.
\emph{Biometrika} \bold{63}, 493--500.

Cox, T.F. (1979)
A method for mapping the dense and sparse regions of
a forest stand.
\emph{Applied Statistics} \bold{28}, 14--19.

Cox, T.F. and Lewis, T. (1976)
A conditioned distance ratio method for analysing spatial patterns.
\emph{Biometrika} \bold{63}, 483--492.

Diggle, P.J. (1979a)
The detection of random heterogeneity in plant populations.
\emph{Biometrics} \bold{33}, 390--394.

Diggle, P.J. (1979b)
Statistical methods for spatial point patterns in ecology.
\emph{Spatial and temporal analysis in ecology}.
R.M. Cormack and J.K. Ord (eds.)
Fairland: International Co-operative Publishing House.
pages 95--150.

Diggle, P.J. (1981)
Some graphical methods in the analysis of spatial point patterns.
In \emph{Interpreting Multivariate Data}.
V. Barnett (eds.) John Wiley and Sons. Pages 55--73.

Diggle, P.J. (1983)
\emph{Statistical analysis of spatial point patterns}.
Academic Press.

Gerrard, D.J. (1969)
Competition quotient: a new measure of the competition
affecting individual forest trees.
Research Bulletin 20, Agricultural Experiment Station,
Michigan State University.

Lotwick, H.W. (1981)
\emph{Spatial stochastic point processes}.
PhD thesis, University of Bath, UK.

Ord, J.K. (1978) 
How many trees in a forest?
\emph{Mathematical Scientist} \bold{3}, 23--33.

}
\keyword{datasets}
\keyword{spatial}
