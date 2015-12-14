\name{Tstat}
\alias{Tstat}
\title{
  Third order summary statistic
}
\description{
  Computes the third order summary statistic \eqn{T(r)}
  of a spatial point pattern.
}
\usage{
Tstat(X, ..., r = NULL, rmax = NULL,
    correction = c("border", "translate"), ratio = FALSE, verbose=TRUE)
}
\arguments{
  \item{X}{The observed point pattern, 
    from which an estimate of \eqn{T(r)} will be computed.
    An object of class \code{"ppp"}, or data
    in any format acceptable to \code{\link{as.ppp}()}.
  }
  \item{\dots}{Ignored.}
  \item{r}{
    Optional. Vector of values for the argument \eqn{r} at which \eqn{T(r)} 
    should be evaluated. Users are advised \emph{not} to specify this
    argument; there is a sensible default.
  }
  \item{rmax}{
    Optional. Numeric. The maximum value of \eqn{r} for which
    \eqn{T(r)} should be estimated.
  }
  \item{correction}{
    Optional. A character vector containing any selection of the
    options \code{"none"}, \code{"border"}, \code{"bord.modif"},
    \code{"translate"}, \code{"translation"}, or \code{"best"}.
    It specifies the edge correction(s) to be applied.
  }
  \item{ratio}{
    Logical. 
    If \code{TRUE}, the numerator and denominator of
    each edge-corrected estimate will also be saved,
    for use in analysing replicated point patterns.
  }
  \item{verbose}{
    Logical. If \code{TRUE}, an estimate of the computation time
    is printed.
  }
}
\details{
  This command calculates the 
  third-order summary statistic \eqn{T(r)} for a spatial point patterns,
  defined by Schladitz and Baddeley (2000).

  The definition of \eqn{T(r)} is similar to the definition of Ripley's
  \eqn{K} function \eqn{K(r)}, except that \eqn{K(r)} counts pairs of
  points while \eqn{T(r)} counts triples of points. 
  Essentially \eqn{T(r)} is a rescaled cumulative
  distribution function of the diameters of triangles in the
  point pattern. The diameter of a triangle is the length of its
  longest side.
}
\section{Computation time}{
  If the number of points is large, the algorithm can take a very long time
  to inspect all possible triangles. A rough estimate
  of the total computation time will be printed at the beginning
  of the calculation. If this estimate seems very large,
  stop the calculation using the user interrupt signal, and
  call \code{Tstat} again, using \code{rmax} to restrict the
  range of \code{r} values,
  thus reducing the number of triangles to be inspected.
}
\value{
  An object of class \code{"fv"}, see \code{\link{fv.object}},
  which can be plotted directly using \code{\link{plot.fv}}.
}
\references{
  Schladitz, K. and Baddeley, A. (2000)
  A third order point process characteristic.
  \emph{Scandinavian Journal of Statistics} \bold{27} (2000) 657--671.
}
\seealso{
  \code{\link{Kest}}
}
\examples{
  plot(Tstat(redwood))
}
\author{Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
}
\keyword{spatial}
\keyword{nonparametric}
