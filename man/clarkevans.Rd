\name{clarkevans}
\alias{clarkevans}
\title{Clark and Evans Aggregation Index}
\description{
  Computes the Clark and Evans aggregation index
  \eqn{R} for a spatial point pattern.
}
\usage{
clarkevans(X, correction=c("none", "Donnelly", "cdf"),
              clipregion=NULL)
}
\arguments{
  \item{X}{
    A spatial point pattern (object of class \code{"ppp"}).
  }
  \item{correction}{
    Character vector.
    The type of edge correction(s) to be applied.
  }
  \item{clipregion}{
    Clipping region for the guard area correction.
    A window (object of class \code{"owin"}).
    See Details.
  }
}
\details{
  The Clark and Evans (1954) aggregation index \eqn{R} is a crude
  measure of clustering or ordering of a point pattern.
  It is the ratio of the observed mean nearest neighbour distance
  in the pattern to that expected for a Poisson point process
  of the same intensity.
  A value \eqn{R>1} suggests ordering, while \eqn{R<1} suggests
  clustering.

  Without correction for edge effects, the value of \code{R} will be
  positively biased. Edge effects arise because, for a point of \code{X}
  close to the edge of the window, the true nearest neighbour may
  actually lie outside the window. Hence observed nearest neighbour
  distances tend to be larger than the true nearest neighbour distances.

  The argument \code{correction} specifies an edge correction
  or several edge corrections to be applied. It is a character vector
  containing one or more of the options
  \code{"none"}, \code{"Donnelly"}, \code{"guard"} and \code{"cdf"}
  (which are recognised by partial matching).
  These edge corrections are:
  \describe{
    \item{"none":}{
      No edge correction is applied. 
    }
    \item{"Donnelly":}{
      Edge correction of Donnelly (1978), available for rectangular
      windows only.
      The theoretical expected value of mean nearest neighbour distance
      under a Poisson process is adjusted for edge effects by the 
      edge correction of Donnelly (1978). The value of \eqn{R} is the
      ratio of the observed mean nearest neighbour distance to this
      adjusted theoretical mean.
    }
    \item{"guard":}{
      Guard region or buffer area method.
      The observed mean nearest neighbour distance
      for the point pattern \code{X}
      is re-defined by averaging only over those points of \code{X}
      that fall inside the sub-window \code{clipregion}.
    }
    \item{"cdf":}{
      Cumulative Distribution Function method.
      The nearest neighbour distance distribution function \eqn{G(r)}
      of the stationary point process is estimated by \code{\link{Gest}}
      using the Kaplan-Meier type edge correction. Then the mean of the
      distribution is calculated from the cdf. 
    }
  }
  Alternatively \code{correction="all"} selects all options.

  If the argument \code{clipregion} is given, then the selected
  edge corrections will be assumed to include \code{correction="guard"}.

  To perform a test based on the Clark-Evans index,
  see \code{\link{clarkevans.test}}.
}
\value{
  A numeric value, or a numeric vector with named components
  \item{naive}{\eqn{R} without edge correction}
  \item{Donnelly}{\eqn{R} using Donnelly edge correction}
  \item{guard}{\eqn{R} using guard region}
  \item{cdf}{\eqn{R} using cdf method}
  (as selected by \code{correction}). The value of the \code{Donnelly}
  component will be \code{NA} if the window of \code{X} is not a rectangle.
}
\references{
  Clark, P.J. and Evans, F.C. (1954)
  Distance to nearest neighbour as a measure of spatial
  relationships in populations \emph{Ecology} \bold{35},
  445--453.

  Donnelly, K. (1978) Simulations to determine the variance
  and edge-effect of total nearest neighbour distance.
  In I. Hodder (ed.) \emph{Simulation studies in archaeology},
  Cambridge/New York: Cambridge University Press, pp 91--95.
}
\author{
John Rudge
\email{rudge@esc.cam.ac.uk}
with modifications by
  \adrian
  
  
}
\seealso{
  \code{\link{clarkevans.test}},
  \code{\link{hopskel}},
  \code{\link{nndist}},
  \code{\link{Gest}}
}
\examples{
  # Example of a clustered pattern
  clarkevans(redwood)

  # Example of an ordered pattern
  clarkevans(cells)

  # Random pattern
  X <- rpoispp(100)
  clarkevans(X)

  # How to specify a clipping region
  clip1 <- owin(c(0.1,0.9),c(0.1,0.9))
  clip2 <- erosion(Window(cells), 0.1)
  clarkevans(cells, clipregion=clip1)
  clarkevans(cells, clipregion=clip2)
}
\keyword{spatial}
\keyword{nonparametric}
