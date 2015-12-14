\name{fryplot}
\alias{fryplot}
\alias{frypoints}
\title{Fry Plot of Point Pattern}
\description{
  Displays the Fry plot (Patterson plot) of a spatial point pattern. 
}
\usage{
fryplot(X, ..., width=NULL, from=NULL, to=NULL, axes=FALSE)
frypoints(X, from=NULL, to=NULL, dmax=Inf)
}
\arguments{
  \item{X}{
    A point pattern (object of class \code{"ppp"}) or something
    acceptable to \code{\link{as.ppp}}.
  }
  \item{\dots}{Optional arguments to control the appearance of the
    plot.}
  \item{width}{Optional parameter indicating the width of a box
    for a zoomed-in view of the Fry plot near the origin.}
  \item{from,to}{
    Optional. Subset indices specifying which points of \code{X}
    will be considered when forming the vectors (drawn from each point of
    \code{from}, to each point of \code{to}.)
  }
  \item{axes}{
    Logical value indicating whether to draw axes,
    crossing at the origin.
  }
  \item{dmax}{
    Maximum distance between points. Pairs at greater distances
    do not contribute to the result. The default means there is
    no maximum distance.
  }
}
\details{
  The function \code{fryplot} generates a Fry plot (or Patterson plot);
  \code{frypoints} returns the points of the Fry plot as a point pattern
  dataset.

  Fry (1979) and Hanna and Fry (1979) introduced a manual graphical method for
  investigating features of a spatial point pattern of mineral deposits.
  A transparent sheet, marked
  with an origin or centre point, is placed over the point pattern.
  The transparent sheet is shifted so that the origin lies over one of
  the data points, and the positions of all the \emph{other} data points
  are copied onto the transparent sheet. This procedure is repeated for
  each data point in turn. The resulting plot (the Fry plot)
  is a pattern of \eqn{n(n-1)} points, where \eqn{n} is the original number
  of data points. This procedure was previously proposed by
  Patterson (1934, 1935) for studying inter-atomic distances in
  crystals, and is also known as a Patterson plot.

  The function \code{fryplot} generates the Fry/Patterson plot.
  Standard graphical parameters
  such as \code{main}, \code{pch},
  \code{lwd}, \code{col}, \code{bg}, \code{cex} can be used to control
  the appearance of the plot.
  To zoom in (to view only a subset of the Fry plot at higher
  magnification), use the argument \code{width} to specify the width
  of a rectangular field of view centred at the origin, or the standard
  graphical arguments \code{xlim} and \code{ylim} to specify another
  rectangular field of view. (The actual field of view may be slightly
  larger, depending on the graphics device.)

  The function \code{frypoints} returns the points of the Fry
  plot as a point pattern object. There may be a large number of points
  in this pattern, so this function should be used only if further
  analysis of the Fry plot is required.

  Fry plots are particularly useful for recognising anisotropy in
  regular point patterns. A void around the origin in the Fry plot
  suggests regularity (inhibition between points) and the shape of the
  void gives a clue to anisotropy in the pattern. Fry plots are also
  useful for detecting periodicity or rounding of the spatial
  coordinates. 

  In mathematical terms, the Fry plot of a point pattern \code{X}
  is simply a plot of the vectors \code{X[i] - X[j]} connecting all
  pairs of distinct points in \code{X}. 

  The Fry plot is related to the \eqn{K} function (see
  \code{\link{Kest}}) and the reduced second moment measure
  (see \code{\link{Kmeasure}}). For example, the number 
  of points in the Fry plot lying within a circle of given radius
  is an unnormalised and uncorrected version of the \eqn{K} function.
  The Fry plot has a similar appearance to the plot of the
  reduced second moment measure \code{\link{Kmeasure}} when the
  smoothing parameter \code{sigma} is very small.
  
  The Fry plot does not adjust for the effect
  of the size and shape of the sampling window. 
  The density of points in the Fry plot tapers off near the edges of the
  plot. This is an edge effect, a consequence of the bounded sampling
  window. In geological applications this is usually not
  important, because interest is focused on the behaviour near the
  origin where edge effects can be ignored.
  To correct for the edge effect, use \code{\link{Kmeasure}} or
  \code{\link{Kest}} or its relatives.
}
\value{
  \code{fryplot} returns \code{NULL}.
  \code{frypoints} returns a point pattern (object of class \code{"ppp"}).
}
\references{
  Fry, N. (1979) 
  Random point distributions and strain measurement in rocks.
  \emph{Tectonophysics} \bold{60}, 89--105.

  Hanna, S.S. and Fry, N. (1979)
  A comparison of methods of strain determination in rocks from
  southwest Dyfed (Pembrokeshire) and adjacent areas.
  \emph{Journal of Structural Geology} \bold{1}, 155--162.

  Patterson, A.L. (1934) A Fourier series method for the determination
  of the component of inter-atomic distances in crystals.
  \emph{Physics Reviews} \bold{46}, 372--376.

  Patterson, A.L. (1935) A direct method for the determination of the
  components of inter-atomic distances in crystals.
  \emph{Zeitschrift fuer Krystallographie} \bold{90}, 517--554.
}
\seealso{
  \code{\link{Kmeasure}},
  \code{\link{Kest}}
}
\examples{
## unmarked data
fryplot(cells)
Y <- frypoints(cells)

## numerical marks
fryplot(longleaf, width=4, axes=TRUE)

## multitype points
fryplot(amacrine, width=0.2,
                  from=(marks(amacrine) == "on"),
                  chars=c(3,16), cols=2:3,
                  main="Fry plot centred at an On-cell")
points(0,0)
}
\author{Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  and Rolf Turner \email{r.turner@auckland.ac.nz}
  
}
\keyword{spatial}
\keyword{nonparametric}
