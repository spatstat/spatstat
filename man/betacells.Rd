\name{betacells}
\alias{betacells}
\docType{data}
\title{Beta Ganglion Cells in Cat Retina} 
\description{
  Point pattern of cells in the retina, each cell classified as `on' or
  `off' and labelled with the cell profile area.
} 
\format{
  \code{betacells} is an object of class \code{"ppp"}
  representing the point pattern of cell locations.
  Entries include
  \tabular{ll}{
    \code{x} \tab Cartesian \eqn{x}-coordinate of cell \cr
    \code{y} \tab Cartesian \eqn{y}-coordinate of cell \cr
    \code{marks} \tab data frame of marks
  }
  Cartesian coordinates are given in microns.

  The data frame of marks has two columns:
  \tabular{ll}{
    \code{type} \tab factor with levels \code{off} and \code{on} \cr
                \tab indicating ``off'' and ``on'' cells\cr
    \code{area} \tab numeric vector giving the \cr
		\tab areas of cell profiles (in square microns)
  }
  See \code{\link{ppp.object}} for details of the format.
}
\usage{data(betacells)}
\source{
  W\"assle et al (1981), Figure 6(a),
  scanned and processed by Stephen Eglen
  \email{S.J.Eglen@damtp.cam.ac.uk}
}
\section{Notes}{
  This is a new, corrected version of the old dataset
  \code{\link{ganglia}}. See below.
  
  These data represent a pattern of beta-type ganglion cells
  in the retina of a cat recorded by W\"assle et al. (1981).
  Beta cells are associated
  with the resolution of fine detail in the cat's visual system.
  They can be classified anatomically as ``on'' or ``off''.
 
  Statistical independence of the arrangement of the
  ``on''- and ``off''-components
  would strengthen the evidence for Hering's (1878) `opponent theory'
  that there are two separate channels for sensing
  ``brightness'' and ``darkness''.
  See W\"assle et al (1981). There is considerable current interest
  in the arrangement of cell mosaics in the retina, see
  Rockhill et al (2000).

  The dataset is a marked point pattern giving the locations,
  types (``on'' or ``off''), and profile areas of beta cells observed
  in a rectangle of dimensions \eqn{750 \times 990}{750 x 990} microns.
  Coordinates are given in microns (thousandths of a millimetre)
  and areas are given in square microns.
  
  The original source is Figure 6 of W\"assle et al (1981),
  which is a manual drawing of the beta mosaic observed
  in a microscope field-of-view of a whole mount of the retina.
  Thus, all beta cells in the retina were effectively projected onto the same
  two-dimensional plane.
  
  The data were scanned in 2004 by Stephen Eglen from
  Figure 6(a) of W\"assle et al (1981).
  Image analysis software was used to identify the soma (cell
  body). The \eqn{x,y} location of each cell was taken to be the
  centroid of the soma. The type of each cell (``on'' or `off'')
  was identified by referring to Figures 6(b) and 6(d).
  The area of each soma (in square microns) was also computed.

  Note that this is a corrected version of the \code{\link{ganglia}}  
  dataset provided in earlier versions of \pkg{spatstat}.
  The earlier data \code{\link{ganglia}} were not faithful to the scale
  in the original paper and contain some scanning errors.
}
\examples{
   plot(betacells)
   area <- marks(betacells)$area
   plot(betacells \%mark\% sqrt(area/pi), markscale=1)
}
\references{
Hering, E. (1878) Zur Lehre von Lichtsinn. Vienna.

Van Lieshout, M.N.M. and Baddeley, A.J. (1999)
Indices of dependence between types in multivariate point patterns.
\emph{Scandinavian Journal of Statistics} \bold{26}, 511--532.

Rockhill, R.L., Euler, T. and Masland, R.H. (2000)
Spatial order within but not between types of retinal neurons.
\emph{Proc. Nat. Acad. Sci. USA} \bold{97}(5), 2303--2307.

W\"assle, H., Boycott, B. B. & Illing, R.-B. (1981).
Morphology and mosaic of on- and off-beta cells in the cat retina and
some functional considerations.
\emph{Proc. Roy. Soc. London Ser. B} \bold{212}, 177--195.
}
\keyword{datasets}
\keyword{spatial}
