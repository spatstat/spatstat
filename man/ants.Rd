\name{ants}
\alias{ants}
\alias{ants.extra}
\docType{data}
\title{
  Harkness-Isham ants' nests data 
}
\description{
  These data give the spatial locations of nests of two species of
  ants, \emph{Messor wasmanni} and \emph{Cataglyphis bicolor},
  recorded by Professor R.D. Harkness at a site in northern Greece,
  and described in Harkness \& Isham (1983). 
  The full dataset (supplied here) has an irregular polygonal boundary,
  while most analyses have been confined to two rectangular
  subsets of the pattern (also supplied here).

  The harvester ant \emph{M. wasmanni} 
  collects seeds for food and builds a nest composed mainly of
  seed husks. \emph{C. bicolor} is a heat-tolerant desert foraging ant which
  eats dead insects and other arthropods. Interest focuses on whether
  there is evidence in the data for intra-species competition
  between \emph{Messor} nests (i.e. competition for resources) and for
  preferential placement of \emph{Cataglyphis} nests in the vicinity of
  \emph{Messor} nests.   

  The full dataset is displayed in Figure 1 of Harkness \& Isham (1983).
  See \bold{Usage} below to produce a comparable plot.
  It comprises 97 nests (68 Messor and 29 Cataglyphis)
  inside an irregular convex polygonal boundary, together with
  annotations showing a foot track through the region,
  the boundary between field and scrub areas inside the
  region, and indicating the two rectangular subregions
  A and B used in their analysis.
  
  Rectangular subsets of the data were analysed by
  Harkness \& Isham (1983), Isham (1984), Takacs \& Fiksel
  (1986),  S\"arkk\"a (1993, section 5.3),
  H\"ogmander and S\"arkk\"a (1999) and Baddeley \& Turner (2000).
  The full dataset (inside its irregular boundary) was first analysed 
  by Baddeley \& Turner (2005b).
  
  The dataset \code{ants} is the full point pattern
  enclosed by the irregular polygonal boundary.
  The \eqn{x} and \eqn{y} coordinates are eastings (E-W) and northings (N-S)
  scaled so that 1 unit equals 0.5 feet. 
  This is a multitype point pattern object, each point carrying a mark
  indicating the ant species (with levels \code{Cataglyphis}
  and \code{Messor}).

  The dataset \code{ants.extra} is a list of auxiliary
  information:
  \describe{
    \item{\code{A} and \code{B}}{The subsets
      of the pattern within the rectangles A and B
      demarcated in Figure 1 of Harkness \& Isham (1983).
      These are multitype point pattern objects.
    }
    \item{\code{trackNE} and \code{trackSW}}{
      coordinates of two straight lines bounding the foot track.
    }
    \item{\code{fieldscrub}}{The endpoints of a straight line
      separating the regions of `field' and `scrub': 
      scrub to the North and field to the South.
    }
    \item{\code{side}}{
      A \code{function(x,y)} that determines whether the location
      \code{(x,y)} is in the scrub or the field. The function can be applied
      to numeric vectors \code{x} and \code{y}, and returns a factor
      with levels \code{"scrub"} and \code{"field"}.
      This function is useful as a spatial covariate.
    }
    \item{\code{plotit}}{A function which
      produces a plot of the full dataset.
    }
  }
} 
\format{
  \code{ants} is an object of class \code{"ppp"}
  representing the full point pattern of ants' nests.
  See \code{\link{ppp.object}} for details of the format.
  The coordinates are scaled so that 1 unit equals 0.5 feet.
  The points are marked by species (with levels \code{Cataglyphis}
  and \code{Messor}).

  \code{ants.extra} is a list with entries
  \describe{
    \item{A}{point pattern of class \code{"ppp"}}
    \item{B}{point pattern of class \code{"ppp"}}
    \item{trackNE}{data in format \code{list(x=numeric(2),y=numeric(2))}
      giving the two endpoints of line markings}
    \item{trackSW}{data in format \code{list(x=numeric(2),y=numeric(2))}
      giving the two endpoints of line markings}
    \item{fieldscrub}{data in format \code{list(x=numeric(2),y=numeric(2))}
      giving the two endpoints of line markings}
    \item{side}{Function with arguments \code{x,y}}
    \item{plotit}{Function}
  }
}
\usage{data(ants)}
\examples{

  # Equivalent to Figure 1 of Harkness and Isham (1983)

  data(ants)
  ants.extra$plotit()

  # Data in subrectangle A, rotated 
  # Approximate data used by Sarkka (1993)

  angle <- atan(diff(ants.extra$fieldscrub$y)/diff(ants.extra$fieldscrub$x))
  plot(rotate(ants.extra$A, -angle))

  # Approximate window used by Takacs and Fiksel (1986)

  tfwindow <- bounding.box(ants$window)
  antsTF <- ppp(ants$x, ants$y, window=tfwindow)
  plot(antsTF)
}
\source{Harkness and Isham (1983).
  Nest coordinates kindly provided by Prof Valerie Isham.
  Polygon coordinates digitised by Adrian Baddeley from
  a reprint of Harkness \& Isham (1983).
}
\references{
Baddeley, A. and Turner, R. (2000)
Practical maximum pseudolikelihood for spatial point patterns.
\emph{Australian and New Zealand Journal of Statistics}
\bold{42}, 283--322.
 
  Baddeley, A. and Turner, R. (2005a)
  Spatstat: an R package for analyzing spatial point patterns.
  \emph{Journal of Statistical Software} \bold{12}:6, 1--42.
  URL: \code{www.jstatsoft.org}, ISSN: 1548-7660.

  Baddeley, A. and Turner, R. (2005b)
  Modelling spatial point patterns in R.
  In: A. Baddeley, P. Gregori, J. Mateu, R. Stoica, and D. Stoyan,
  editors, \emph{Case Studies in Spatial Point Pattern Modelling},
  Lecture Notes in Statistics number 185. Pages 23--74.
  Springer-Verlag, New York, 2006. 
  ISBN: 0-387-28311-0.  

Harkness, R.D. and Isham, V. (1983)
A bivariate spatial point pattern of ants' nests.
\emph{Applied Statistics} \bold{32}, 293--303.

\ifelse{latex}{\out{H\"ogmander}}{Hogmander}, H. and
\ifelse{latex}{\out{S\"arkk\"a}}{Sarkka}, A. (1999)
Multitype spatial point patterns with hierarchical interactions.
\emph{Biometrics} \bold{55}, 1051--1058.

Isham, V.S. (1984)
Multitype Markov point processes: some approximations.
\emph{Proceedings of the Royal Society of London, Series A},
\bold{391}, 39--53.

Takacs, R. and Fiksel, T. (1986)
Interaction pair-potentials for a system of ants' nests.
\emph{Biometrical Journal} \bold{28}, 1007--1013.

\ifelse{latex}{\out{S\"arkk\"a}}{Sarkka}, A. (1993)
\emph{Pseudo-likelihood approach for pair potential
  estimation of Gibbs processes}.
Number 22 in \ifelse{latex}{\out{Jyv\"askyl\"a}}{Jyvaskyla}
Studies in Computer Science, Economics and Statistics.
University of \ifelse{latex}{\out{Jyv\"askyl\"a}}{Jyvaskyla}, Finland.
}
\keyword{datasets}
\keyword{spatial}
