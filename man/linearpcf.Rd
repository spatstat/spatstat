\name{linearpcf}
\alias{linearpcf}
\title{
  Linear Pair Correlation Function
}
\description{
  Computes an estimate of the linear pair correlation function
  for a point pattern on a linear network.
}
\usage{
linearpcf(X, r=NULL, ..., correction="Ang")
}
\arguments{
  \item{X}{
    Point pattern on linear network (object of class \code{"lpp"}).
  }
  \item{r}{
    Optional. Numeric vector of values of the function argument \eqn{r}.
    There is a sensible default.
  }
  \item{\dots}{
    Arguments passed to \code{\link{density.default}}
    to control the smoothing.
  }
  \item{correction}{
    Geometry correction.
    Either \code{"none"} or \code{"Ang"}. See Details.
  }
}
\details{
  This command computes the linear pair correlation function
  from point pattern data on a linear network.

  The pair correlation function is estimated from the
  shortest-path distances between each pair of data points,
  using the fixed-bandwidth kernel smoother
  \code{\link{density.default}}, 
  with a bias correction at each end of the interval of \eqn{r} values.
  To switch off the bias correction, set \code{endcorrect=FALSE}.
  
  If \code{correction="none"}, the calculations do not include
  any correction for the geometry of the linear network. The result is
  an estimate of the first derivative of the 
  network \eqn{K} function defined by Okabe and Yamada (2001).

  If \code{correction="Ang"}, the pair counts are weighted using
  Ang's correction (Ang, 2010). The result is an estimate of the
  pair correlation function in the linear network.

}
\value{
  Function value table (object of class \code{"fv"}).
}
\author{
  Ang Qi Wei \email{aqw07398@hotmail.com} and
  Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
}
\references{
  Ang, Q.W. (2010) Statistical methodology for spatial point patterns
  on a linear network. MSc thesis, University of Western Australia.
  
  Ang, Q.W., Baddeley, A. and Nair, G. (2012)
  Geometrically corrected second-order analysis of 
  events on a linear network, with applications to
  ecology and criminology.
  \emph{Scandinavian Journal of Statistics} \bold{39}, 591--617.
  
  Okabe, A. and Yamada, I. (2001) The K-function method on a network and
  its computational implementation. \emph{Geographical Analysis}
  \bold{33}, 271-290.
}
\seealso{
  \code{\link{linearK}},
  \code{\link{linearpcfinhom}},
  \code{\link{lpp}}
}
\examples{
  data(simplenet)
  X <- rpoislpp(5, simplenet)
  linearpcf(X)
  linearpcf(X, correction="none")
}
\keyword{spatial}
\keyword{nonparametric}
