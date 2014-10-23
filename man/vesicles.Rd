\name{vesicles}
\alias{vesicles}
\alias{vesicles.extra}
\docType{data}
\title{
  Vesicles Data
}
\description{
  Point pattern of synaptic vesicles observed in rat brain tissue.
}
\usage{data(vesicles)}
\format{
  The dataset \code{vesicles} is a point pattern
  (object of class \code{"ppp"}) representing the location
  of the synaptic vesicles. The window of the point pattern
  represents the region of presynapse where synaptic vesicles were
  observed in this study.
  There is a hole in the window, representing the region occupied by
  mitochondria, where synaptic vesicles do not occur.

  The dataset \code{vesicles.extra} is a list with entries
  \tabular{ll}{
    \code{presynapse}\tab outer polygonal boundary of presynapse \cr
    \code{mitochondria} \tab polygonal boundary of mitochondria \cr
    \code{mask} \tab binary mask representation of vesicles window \cr
    \code{activezone} \tab line segment pattern representing \cr
                      \tab the active zone.
  }
}
\details{
  As part of a study on the effects of stress on brain function,
  Khanmohammadi et al (2014) analysed the spatial pattern of
  synaptic vesicles in 45-nanometre-thick sections of
  rat brain tissue visualised in transmission electron microscopy.
  
  To investigate the influence of stress, Khanmohammadi et al (2014)
  study the distribution of the synaptic vesicles in
  the pre-synaptic neuron in relation to the active zone of the
  presynaptic membrane. They hypothesize that the synaptic vesicle density
  is a decreasing function of distance to the active zone. 

  The boundaries for the active zone, mitochondria, pre- and post
  synaptic terminals, and the centre of the synaptic vesicles
  were annotated by hand on the image.
}
\section{Raw Data}{
  For demonstration and training purposes,
  the raw data files for this dataset are also
  provided in the \pkg{spatstat} package installation:
  \tabular{ll}{
    \code{vesicles.txt}\tab spatial locations of vesicles \cr
    \code{presynapse.txt}\tab vertices of \code{presynapse} \cr
    \code{mitochondria.txt} \tab vertices of \code{mitochondria} \cr
    \code{vesiclesimage.tif} \tab greyscale microscope image \cr
    \code{vesiclesmask.tif} \tab binary image of \code{mask} \cr
    \code{activezone.txt} \tab coordinates of \code{activezone} 
  }
  The files are in the folder \code{rawdata/vesicles} in the
  \pkg{spatstat} installation directory. The precise location of the
  files can be obtained using \code{\link[base]{system.file}}, as shown
  in the examples.  
}
\source{
  Nicoletta Nava, Mahdieh Khanmohammadi and Jens Randel Nyengaard.
  Experiment performed by Nicoletta Nava at the
  Stereology and Electron Microscopy Laboratory, Aarhus University,
  Denmark. Images were annotated by Mahdieh Khanmohammadi
  at the Department of Computer Science, University of Copenhagen.
  Jens Randel Nyengaard provided supervision and guidance,
  and curated the data.
}
\references{
  Khanmohammadi, M., Waagepetersen, R., Nava, N.,
  Nyengaard, J.R. and Sporring, J. (2014)
  Analysing the distribution of synaptic vesicles using a
  spatial point process model.
  \emph{5th ACM Conference on Bioinformatics, Computational Biology
  and Health Informatics}, Newport Beach, CA, USA, September 2014.
}
\examples{
plot(vesicles)
with(vesicles.extra,
     plot(activezone, add=TRUE, col="red"))

## read coordinates of vesicles from raw data, for training purposes
vf <- system.file("rawdata/vesicles/vesicles.txt", package="spatstat")
vdf <- read.table(vf, header=TRUE)
}
\keyword{spatial}
\keyword{datasets}
