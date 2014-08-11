\name{gorillas}
\alias{gorillas}
\alias{gorillas.extra}
\docType{data}
\title{
  Gorilla Nesting Sites
}
\description{
  Locations of nesting sites of gorillas, and associated covariates,
  in a National Park in Cameroon.
}
\usage{data(gorillas)}
\format{
  \code{gorillas} is a marked point pattern (object
  of class \code{"ppp"}) representing nest site locations.

  \code{gorillas.extra} is a named list of 7 pixel images (objects of
  class \code{"im"}) containing spatial covariates.
  It also belongs to the class \code{"listof"}.

  All spatial coordinates are in metres.
  The coordinate reference system is \code{WGS_84_UTM_Zone_32N}.
}
\details{
  These data come from a study of gorillas in the Kagwene Gorilla Sanctuary,
  Cameroon, by the Wildlife Conservation Society
  Takamanda-Mone Landscape Project (WCS-TMLP). A detailed description
  and analysis of the data is reported in Funwi-Gabga and Mateu (2012).

  The dataset \code{gorillas} is a marked point pattern
  (object of class \code{"ppp"})
  giving the spatial locations of 647 nesting sites of gorilla groups
  observed in the sanctuary over time.
  Locations are given as UTM (Zone 32N) coordinates in metres.
  The observation window is the boundary of the sanctuary, represented
  as a polygon. Marks attached to the points are:
  
  \describe{
    \item{group}{Identifier of the gorilla group
      that constructed the nest site:
      a categorical variable with values \code{major} or \code{minor}.
    }
    \item{season}{Season in which data were collected:
      categorical, either \code{rainy} or \code{dry}.
    }
    \item{date}{
      Day of observation. A value of class \code{"Date"}.
    }
  }

  The accompanying dataset \code{gorillas.extra} contains
  spatial covariate information. It is a named list containing
  seven pixel images (objects of class \code{"im"}) giving the values of
  seven covariates over the study region. It also belongs
  to the class \code{"listof"} so that it can be plotted.
  The component images are:
  
  \describe{
    \item{aspect}{
      Compass direction of the terrain slope.
      Categorical, with levels
      \code{N}, 
      \code{NE}, 
      \code{E}, 
      \code{SE}, 
      \code{S}, 
      \code{SW}, 
      \code{W} and 
      \code{NW}.
    }
    \item{elevation}{
      Digital elevation of terrain, in metres.
    }
    \item{heat}{
      Heat Load Index at each point on the surface (Beer's aspect),
      discretised. Categorical with values \code{Warmest}
      (Beer's aspect between 0 and 0.999),
      \code{Moderate} (Beer's aspect between 1 and 1.999),
      \code{Coolest} (Beer's aspect equals 2).
    }
    \item{slopeangle}{
      Terrain slope, in degrees.
    }
    \item{slopetype}{
      Type of slope.
      Categorical, with values
      \code{Valley}, \code{Toe} (toe slope), \code{Flat},
      \code{Midslope}, \code{Upper} and \code{Ridge}.
    }
    \item{vegetation}{
      Vegetation or cover type.
      Categorical, with values
      \code{Disturbed} (highly disturbed forest), \code{Colonising}
      (colonising forest), \code{Grassland} (savannah),
      \code{Primary} (primary forest), \code{Secondary} (secondary forest), and
      \code{Transition} (transitional vegetation).
    }
    \item{waterdist}{
      Euclidean distance from nearest water body, in metres.
    }
  }

  For further information see Funwi-Gabga and Mateu (2012).
}
\source{
  Field data collector: Wildlife Conservation Society Takamanda-Mone
  Landscape Project (WCS-TMLP).
  \emph{Please acknowledge WCS-TMLP in any use of these data.}

  Data kindly provided by 
  Funwi-Gabga Neba, Data Coordinator of A.P.E.S.
  Database Project, Department of Primatology,
  Max Planck Institute for Evolutionary Anthropology, Leipzig, Germany.

  The collaboration of Prof Jorge Mateu, Universitat Jaume I, Castellon, Spain
  is gratefully acknowledged.
}
\references{
  Funwi-Gabga, N. (2008)
  \emph{A pastoralist survey and fire impact assessment
  in the Kagwene Gorilla Sanctuary, Cameroon}. M.Sc. thesis,
  Geology and Environmental Science, University of Buea.

  Funwi-Gabga, N. and Mateu, J. (2012)
  Understanding the nesting spatial behaviour of gorillas
  in the Kagwene Sanctuary, Cameroon.
  \emph{Stochastic Environmental Research and Risk Assessment},
  in press.
}
\examples{
  summary(gorillas)
  plot(gorillas)
  plot(gorillas.extra)
}
\keyword{datasets}
