\name{ganglia}
\alias{ganglia}
\docType{data}
\title{Beta Ganglion Cells in Cat Retina, Old Version} 
\description{
  Point pattern of retinal ganglion cells identified as `on' or `off'.
  A marked point pattern.
} 
\format{
  An object of class \code{"ppp"}
  representing the point pattern of cell locations.
  Entries include
  \tabular{ll}{
    \code{x} \tab Cartesian \eqn{x}-coordinate of cell \cr
    \code{y} \tab Cartesian \eqn{y}-coordinate of cell \cr
    \code{marks} \tab factor with levels \code{off} and \code{on} \cr
                 \tab indicating ``off'' and ``on'' cells
  }
  See \code{\link{ppp.object}} for details of the format.
}
\usage{data(ganglia)}
\source{W\"assle et al (1981), data supplied by Marie-Colette van
  Lieshout and attributed to Peter Diggle}
\section{Notes}{
  \bold{Important: these data are INCORRECT.  See below.}
  
  The data represent a pattern of beta-type ganglion cells in the
  retina of a cat recorded in Figure 6(a) of W\"assle et al. (1981).

  The pattern was first analysed by W\"assle et al (1981) using
  nearest neighbour distances. The data used in their analysis
  are not available.

  The present dataset \code{\link{ganglia}} was
  scanned from Figure 6(a) of W\"assle et al (1981)
  in the early 1990's, but we have no further information.
  This dataset is the one analysed by Van Lieshout and Baddeley (1999)
  using multitype J functions, and by Stoyan (1995) using second
  order methods (pair correlation and mark correlation).
  
  It has now been discovered that these data are \bold{incorrect}.
  They are not faithful to the scale in Figure 6 of W\"assle et al (1981),
  and they contain some scanning errors.
  Hence they should not be used to address the original scientific question.
  They have been retained only for comparison with other analyses
  in the statistical literature.

  A new, corrected dataset, scanned from the original microscope image,
  has been provided under the name \code{\link{betacells}}.  Use that
  dataset for any further study.
}
\section{Warnings}{
  These data are incorrect.
  Use the new corrected dataset \code{\link{betacells}}.  
}
\references{
  Stoyan, D. (1995) Personal communication.
  
  Van Lieshout, M.N.M. and Baddeley, A.J. (1999) Indices of
  dependence between types in multivariate point patterns.
  \emph{Scandinavian Journal of Statistics} \bold{26}, 511--532.

  W\"assle, H., Boycott, B. B. & Illing, R.-B. (1981).  Morphology
  and mosaic of on- and off-beta cells in the cat retina and some
  functional considerations.  \emph{Proc. Roy. Soc. London Ser. B}
  \bold{212}, 177--195.
}
\keyword{datasets}
\keyword{spatial}


