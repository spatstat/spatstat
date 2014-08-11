\name{clmfires}
\alias{clmfires}
\alias{clmcov100}
\alias{clmcov200}
\docType{data}
\title{
  Castilla-La Mancha Forest Fires
}
\description{
  This dataset is a record of forest fires in the Castilla-La Mancha
  region of Spain between 1998 and 2007.  This region is approximately
  400 by 400 kilometres.  The coordinates are recorded in kilometres.

  The dataset \code{clmfires} is a point pattern (object of class
  \code{"ppp"}) containing the spatial coordinates of each fire,
  with marks containing information about each fire.  There are 4
  columns of marks:
  \tabular{ll}{
    \code{cause} \tab cause of fire (see below) \cr
    \code{burnt.area} \tab total area burned, in hectares \cr
    \code{date} \tab the date of fire as an object of class \code{Date} \cr
    \code{julian.date} \tab number of days elapsed since 1 January 1998 \cr
  }
  The \code{cause} of the fire is a factor with the levels
  \code{lightning}, \code{accident} (for accidents or negligence),
  \code{intentional} (for intentionally started fires) and
  \code{other} (for other causes including unknown cause).

  The format of \code{date} is \dQuote{Year-month-day}, e.g.
  \dQuote{2005-07-14} (i.e. 14 July, 2005).

  Accompanying this point pattern, there are two datasets
  \code{clmcov100} and \code{clmcov200} containing covariate
  information for the entire Castilla-La Mancha region. Each
  of these two datasets is a list of four images (objects of
  class \code{"im"}) named \code{elevation}, \code{orientation},
  \code{slope} and \code{landuse}.  The \code{landuse} image is
  factor-valued with the factor having levels \code{urban},
  \code{farm} (for farms or orchards), \code{meadow},
  \code{denseforest} (for dense forest), \code{conifer} (for conifer
  forest or plantation), \code{mixedforest}, \code{grassland},
  \code{bush}, \code{scrub} and \code{artifgreen} for artificial
  greens such as golf courses.

  These images (effectively) provide values for the four
  covariates at every location in the study area. The images in
  \code{clmcov100} are 100 by 100 pixels in size, while those in
  \code{clmcov200} are 200 by 200 pixels.  For easy handling,
  \code{clmcov100} and \code{clmcov200} also belong to the
  class \code{"listof"} so that they can be plotted and printed
  immediately.
} 
\format{
  \code{clmfires} is a marked point pattern (object of class \code{"ppp"}). 
  See \code{\link[spatstat]{ppp.object}}.

  \code{clmcov100} and \code{clmcov200} are lists of pixel images
  (objects of class \code{"im"}).
}
\section{Remark}{
The precision with which the coordinates of the locations of the
fires changed between 2003 and 2004.  From 1998 to 2003 many of
the locations were recorded as the centroid of the corresponding
\dQuote{district unit}; the rest were recorded as exact UTM
coordinates of the centroids of the fires.  In 2004 the system
changed and the exact UTM coordinates of the centroids of the fires
were used for \emph{all} fires.  There is thus a strongly apparent
\dQuote{gridlike} quality to the fire locations for the years 1998
to 2003.

There is however no actual duplication of points in the 1998 to 2003
patterns due to \dQuote{jittering} having been applied in order to
avoid such duplication.  It is not clear just \emph{how} the fire
locations were jittered.  It seems unlikely that the jittering was
done using the \code{jitter()} function from \code{R} or the
\pkg{spatstat} function \code{\link{rjitter}}.

Of course there are many sets of points which are \emph{virtually}
identical, being separated by distances induced by the jittering.
Typically these distances are of the order of 40 metres which
is unlikely to be meaningful on the scale at which forest fires
are observed.

Caution should therefore be exercised in any analyses of the patterns
for the years 1998 to 2003.
}

\usage{data(clmfires)}
\examples{
plot(clmfires, which.marks="cause", cols=2:5, cex=0.25)
plot(clmcov100)
# Split the clmfires pattern by year and plot the first and last years:
yr  <- factor(format(marks(clmfires)$date,format="\%Y"))
X   <- split(clmfires,f=yr)
fAl <- c("1998","2007")
plot(X[fAl],use.marks=FALSE,main.panel=fAl,main="")
}
\source{
  Jorge Mateu.
}
\keyword{datasets}
\keyword{spatial}
