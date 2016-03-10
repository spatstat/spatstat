\name{bramblecanes}
\alias{bramblecanes}
\docType{data}
\title{Hutchings' Bramble Canes data}
\description{
  Data giving the locations and ages of bramble canes in a field.
  A marked point pattern. 
} 
\format{
  An object of class \code{"ppp"} 
  representing the point pattern of plant locations.
  Entries include
  \tabular{ll}{
    \code{x} \tab Cartesian \eqn{x}-coordinate of plant \cr
    \code{y} \tab Cartesian \eqn{y}-coordinate of plant \cr
    \code{marks} \tab factor with levels 0,1, 2 indicating age
  }
  See \code{\link{ppp.object}} for details of the format.
}
\usage{data(bramblecanes)}
\source{Hutchings (1979), data published in Diggle (1983)}
\section{Notes}{
These data record the \eqn{(x,y)} locations and ages of bramble canes in a 
field \eqn{9} metres square, rescaled to the unit square. 
The canes were classified according to age as either newly emergent, 
one or two years old. These are encoded as marks 0, 1 and 2 respectively
in the dataset.

The data were recorded and analysed by Hutchings (1979)
and further analysed by Diggle (1981a, 1981b, 1983),
Diggle and Milne (1983), and Van Lieshout and Baddeley (1999).
All analyses found that the pattern of newly emergent canes
exhibits clustering, which Hutchings attributes to ``vigorous
vegetative reproduction''. 
}
\references{
  Diggle, P. J. (1981a)
  Some graphical methods in the analysis of spatial point patterns.
  In \emph{Interpreting multivariate data}, V. Barnett (Ed.)
  John Wiley and Sons. 
 
  Diggle, P. J. (1981b).
  Statistical analysis of spatial point patterns.
  \emph{N.Z. Statist.} \bold{16}, 22--41.
 
  Diggle, P.J. (1983)
  \emph{Statistical analysis of spatial point patterns}.
  Academic Press.

  Diggle, P. J. and Milne, R. K. (1983)
  Bivariate Cox processes: some models for bivariate spatial point patterns.
  \emph{Journal of the Royal Statistical Soc. Series B} 
  \bold{45}, 11--21.
 
  Hutchings, M. J. (1979)
  Standing crop and pattern in pure stands of Mercurialis
  perennis and Rubus fruticosus in mixed deciduous woodland.
  \emph{Oikos} \bold{31}, 351--357.
 
  Van Lieshout, M.N.M. and Baddeley, A.J. (1999)
  Indices of dependence between types in multivariate point patterns.
  \emph{Scandinavian Journal of Statistics} \bold{26}, 511--532.
}
\keyword{datasets}
\keyword{spatial}
 
