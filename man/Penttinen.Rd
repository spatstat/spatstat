\name{Penttinen}
\alias{Penttinen}
\title{Penttinen Interaction}
\description{
  Creates an instance of the Penttinen pairwise interaction
  point process model, which can then be fitted to point pattern data.
}
\usage{
  Penttinen(r)
}
\arguments{
  \item{r}{circle radius}
}
\value{
  An object of class \code{"interact"}
  describing the interpoint interaction
  structure of a point process. 
}
\details{
  Penttinen (1984, Example 2.1, page 18), citing Cormack (1979),
  described the pairwise interaction point process with interaction factor
  \deqn{
    h(d) = e^{\theta A(d)} = \gamma^{A(d)}
  }{
    h(d) = exp(theta * A(d)) = gamma^(A(d))
  }
  between each pair of points separated by a distance $d$.
  Here \eqn{A(d)} is the area of intersection between two discs
  of radius \eqn{r} separated by a distance \eqn{d}, normalised so that
  \eqn{A(0) = 1}.

  The scale of interaction is controlled by the disc radius \eqn{r}:
  two points interact if they are closer than \eqn{2 r}{2 * r} apart.
  The strength of interaction is controlled by the
  canonical parameter \eqn{\theta}{theta}, which 
  must be less than or equal to zero, or equivalently by the
  parameter \eqn{\gamma = e^\theta}{gamma = exp(theta)},
  which must lie between 0 and 1.

  The potential is inhibitory, i.e.\ this model is only appropriate for
  regular point patterns. 
  For \eqn{\gamma=0}{gamma=0} the model is
  a hard core process with hard core diameter \eqn{2 r}{2 * r}.
  For \eqn{\gamma=1}{gamma=1} the model is a Poisson process.

  The irregular parameter
  \eqn{r} must be given in the call to
  \code{Penttinen}, while the
  regular parameter \eqn{\theta}{theta} will be estimated.

  This model can be considered as a pairwise approximation
  to the area-interaction model \code{\link{AreaInter}}.
}
\seealso{
  \code{\link{ppm}},
  \code{\link{ppm.object}},
  \code{\link{Pairwise}},
  \code{\link{AreaInter}}.
}
\examples{
   fit <- ppm(cells ~ 1, Penttinen(0.07))
   fit
   reach(fit) # interaction range is circle DIAMETER
}
\references{
  Cormack, R.M. (1979)
  Spatial aspects of competition between individuals.
  Pages 151--212 in \emph{Spatial and Temporal Analysis in Ecology},
  eds. R.M. Cormack and J.K. Ord, International Co-operative
  Publishing House, Fairland, MD, USA. 

  Penttinen, A. (1984) 
  \emph{Modelling Interaction in Spatial Point Patterns:
  Parameter Estimation by the Maximum Likelihood Method.}
  \ifelse{latex}{\out{Jyv\"askyl\"a}}{Jyvaskyla}
  Studies in Computer Science, Economics and Statistics \bold{7},
  University of \ifelse{latex}{\out{Jyv\"askyl\"a}}{Jyvaskyla}, Finland.
}
\author{
  \spatstatAuthors
}
\keyword{spatial}
\keyword{models}

