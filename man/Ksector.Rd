\name{Ksector}
\alias{Ksector}
\title{Sector K-function}
\description{
  A directional counterpart of Ripley's \eqn{K} function,
  in which pairs of points are counted only when the
  vector joining the pair happens to
  lie in a particular range of angles.
}
\usage{
Ksector(X, begin = 0, end = 2 * pi, \dots,
        r = NULL, breaks = NULL,
        correction = c("border", "isotropic", "Ripley", "translate"),
        ratio = FALSE)
}
\arguments{
  \item{X}{The observed point pattern, 
    from which an estimate of \eqn{K(r)} will be computed.
    An object of class \code{"ppp"}, or data
    in any format acceptable to \code{\link{as.ppp}()}.
  }
  \item{begin,end}{
    Numeric values giving the range of angles inside which
    points will be counted. Angles are measured in radians
    anti-clockwise from the positive \eqn{x}-axis.
  }
  \item{\dots}{Ignored.}
  \item{r}{
    Optional. Vector of values for the argument \eqn{r} at which \eqn{K(r)} 
    should be evaluated. Users are advised \emph{not} to specify this
    argument; there is a sensible default.
  }
  \item{breaks}{
    This argument is for internal use only.
  }
  \item{correction}{
    Optional. A character vector containing any selection of the
    options \code{"none"}, \code{"border"}, \code{"bord.modif"},
    \code{"isotropic"}, \code{"Ripley"}, \code{"translate"},
    \code{"translation"}, \code{"none"}, \code{"good"} or \code{"best"}.
    It specifies the edge correction(s) to be applied.
  }
  \item{ratio}{
    Logical. 
    If \code{TRUE}, the numerator and denominator of
    each edge-corrected estimate will also be saved,
    for use in analysing replicated point patterns.
  }
}
\details{
  This is a directional counterpart of Ripley's \eqn{K} function
  (see \code{\link{Kest}}) in which, instead of counting all
  pairs of points within a specified distance \eqn{r}, we
  count only the pairs \eqn{(x_i, x_j)}{x[i], x[j]}
  for which the vector \eqn{x_j - x_i}{x[j] - x[i]}
  falls in a particular range of angles.

  This can be used to evaluate evidence for anisotropy
  in the point pattern \code{X}.
}
\value{
  An object of class \code{"fv"} containing the estimated
  function.
}
\seealso{
  \code{\link{Kest}}
}
\examples{
 K <- Ksector(swedishpines, 0, pi/2)
 plot(K)
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  Rolf Turner
  \email{r.turner@auckland.ac.nz}
  and Ege Rubak
  \email{rubak@math.aau.dk}
}
\keyword{spatial}
\keyword{nonparametric}
