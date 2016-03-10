\name{redwoodfull}
\alias{redwoodfull}
\alias{redwoodfull.extra}
\docType{data}
\title{
  California Redwoods Point Pattern (Entire Dataset)
}
\description{
  These data represent the locations of 195 seedlings and saplings 
  of California redwood trees in a square sampling region.  
  They were described and analysed by Strauss (1975).
  This is the ``\bold{full}'' dataset; most writers have
  analysed a subset extracted by Ripley (1977)
  which is available as \code{\link{redwood}}.
 
  Strauss (1975) divided the sampling region into two subregions I and II
  demarcated by a diagonal line.  The spatial pattern
  appears to be slightly regular in region I and strongly clustered in 
  region II. 

  Strauss (1975) writes: \dQuote{It was felt that the seedlings
    would be scattered fairly randomly, except that a number of
    tight clusters would form around some of the redwood tree stumps
    present in the plot. A discontinuity in the soil, very roughly
    demarked by the diagonal line in the figure, was expected to cause
    a difference in clustering behaviour between regions I and II. Moreover,
    almost all the redwood stumps were situated in region II.}

  The dataset \code{redwoodfull} contains the full point pattern
  of 195 trees. 
  The window has been rescaled to the unit square.
  Its physical size is approximately 130 feet across.

  The auxiliary information about the subregions is contained in 
  \code{redwoodfull.extra}, which is a list with entries
  \tabular{ll}{
    \code{rdiag}\tab The coordinates of the diagonal boundary\cr
               \tab between regions I and II \cr
    \code{regionI} \tab Region I as a window object \cr
    \code{regionII} \tab Region II as a window object \cr
    \code{regionR} \tab Ripley's subrectangle (approximate) \cr
    \code{plotit}    \tab Function to plot the full data and auxiliary markings
  }

  Ripley (1977) extracted a subset of these data, containing 62 points,
  lying within a square subregion which overlaps regions I and II.
  He rescaled that subset to the unit square. 
  This subset has been re-analysed many times,
  and is the dataset usually known as
  ``the redwood data'' in the spatial statistics literature.
  The exact dataset used by Ripley is supplied in the \pkg{spatstat}
  library as \code{\link{redwood}}.

  The approximate position of the square chosen by Ripley
  within the \code{redwoodfull} pattern
  is indicated by the window \code{redwoodfull.extra$regionR}.
  There are some minor inconsistencies with
  \code{redwood} since it originates from a different digitisation.
} 
\format{
  The dataset \code{redwoodfull} is an object of class \code{"ppp"}
  representing the point pattern of tree locations.
  See \code{\link{ppp.object}} for details of the format of a
  point pattern object.
  The window has been rescaled to the unit square.
  Its physical size is approximately 128 feet across.

  The dataset \code{redwoodfull.extra} is a list with entries
  \tabular{ll}{
    \code{rdiag}\tab coordinates of endpoints of a line,\cr
               \tab in format \code{list(x=numeric(2),y=numeric(2))} \cr
    \code{regionI} \tab a window object \cr
    \code{regionII} \tab a window object \cr
    \code{regionR} \tab a window object \cr
    \code{plotit}    \tab Function with no arguments
  }
}
\usage{data(redwoodfull)}
\examples{
       data(redwoodfull)
       plot(redwoodfull)
       redwoodfull.extra$plotit()
       # extract the pattern in region II 
       redwoodII <- redwoodfull[, redwoodfull.extra$regionII]
}
\source{Strauss (1975). The plot of the data published by Strauss (1975)
  was scanned and digitised by Sandra Pereira, University of
  Western Australia, 2002.
}
\seealso{
  \code{\link{redwood}}
}
\references{
  Diggle, P.J. (1983)
  \emph{Statistical analysis of spatial point patterns}.
  Academic Press.

  Ripley, B.D. (1977)
  Modelling spatial patterns (with discussion).
  \emph{Journal of the Royal Statistical Society, Series B}
  \bold{39}, 172--212.

  Strauss, D.J. (1975)
  A model for clustering.
  \emph{Biometrika} \bold{62}, 467--475.
}
\keyword{datasets}
\keyword{spatial}
