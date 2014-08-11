\name{bronzefilter}
\alias{bronzefilter}
\docType{data}
\title{Bronze gradient filter data}
\description{
  These data represent a spatially inhomogeneous pattern of
  circular section profiles of particles, observed in a
  longitudinal plane section through a gradient sinter
  filter made from bronze powder, prepared by Ricardo Bernhardt, Dresden. 

  The material was produced by sedimentation of bronze powder with varying 
  grain diameter and subsequent sintering,
  as described in Bernhardt et al. (1997). 

  The data are supplied as a marked point pattern of circle centres marked by
  circle radii.
  The coordinates of the centres and the radii are recorded in mm.
  The field of view is an \eqn{18 \times 7}{18 * 7} mm rectangle.

  The data were first analysed by Hahn et al. (1999).
} 
\format{
  An object of class \code{"ppp"}
  representing the point pattern of cell locations.
  Entries include
  \tabular{ll}{
    \code{x} \tab Cartesian \eqn{x}-coordinate of bronze grain profile centre\cr
    \code{y} \tab Cartesian \eqn{y}-coordinate of bronze grain profile centre\cr
    \code{marks} \tab radius of bronze grain profile
  }
  See \code{\link{ppp.object}} for details of the format.
  All coordinates are recorded in mm.
}
\usage{data(bronzefilter)}
\examples{
  data(bronzefilter)
  plot(bronzefilter, markscale=1)
}
\source{
  R.\ Bernhardt (section image), H.\ Wendrock (coordinate measurement). 
  Adjusted, formatted and communicated by U.\ Hahn.
}
\references{
  Bernhardt, R., Meyer-Olbersleben, F. and Kieback, B. (1997)
  Fundamental investigation on the preparation of gradient structures
  by sedimentation of different powder fractions under gravity.
  \emph{Proc. of the 4th Int. Conf. On Composite Engineering,
    July 6--12 1997, ICCE/4}, Hawaii, Ed. David Hui, 147--148.

  Hahn U., Micheletti, A., Pohlink, R., Stoyan D. and Wendrock, H.(1999) 
  Stereological analysis and modelling of gradient structures. 
  \emph{Journal of Microscopy}, \bold{195}, 113--124.
} 
\keyword{datasets}
\keyword{spatial}
