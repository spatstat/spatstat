\name{StraussHard}
\alias{StraussHard}
\title{The Strauss / Hard Core Point Process Model}
\description{
Creates an instance of the ``Strauss/ hard core'' point process model
which can then be fitted to point pattern data.
}
\usage{
  StraussHard(r, hc)
}
\arguments{
  \item{r}{The interaction radius of the Strauss interaction}
  \item{hc}{The hard core distance}
}
\value{
  An object of class \code{"interact"}
  describing the interpoint interaction
  structure of the ``Strauss/hard core''
  process with Strauss interaction radius \eqn{r}
  and hard core distance \code{hc}.
}
\details{
  A Strauss/hard core process with interaction radius \eqn{r},
  hard core distance \eqn{h < r}, and 
  parameters \eqn{\beta}{beta} and \eqn{\gamma}{gamma},
  is a pairwise interaction point process
  in which
  \itemize{
    \item distinct points are not allowed to come closer
      than a distance \eqn{h} apart
    \item each pair of points closer than \eqn{r} units apart
      contributes a factor \eqn{\gamma}{gamma} to the probability density.
  }
  This is a hybrid of the Strauss process and the hard core process.

  The probability density is zero if any pair of points
  is closer than \eqn{h} units apart, and otherwise equals
  \deqn{
    f(x_1,\ldots,x_n) =
    \alpha \beta^{n(x)} \gamma^{s(x)}
  }{
    f(x_1,\ldots,x_n) =
    alpha . beta^n(x) gamma^s(x)
  }
  where \eqn{x_1,\ldots,x_n}{x[1],\ldots,x[n]} represent the 
  points of the pattern, \eqn{n(x)} is the number of points in the
  pattern, \eqn{s(x)} is the number of distinct unordered pairs of
  points that are closer than \eqn{r} units apart,
  and \eqn{\alpha}{alpha} is the normalising constant.

  The interaction parameter \eqn{\gamma}{gamma} may take any
  positive value (unlike the case for the Strauss process).
  If \eqn{\gamma = 1}{gamma = 1}, the process reduces to a classical
  hard core process.
  If \eqn{\gamma < 1}{gamma < 1},
  the model describes an ``ordered'' or ``inhibitive'' pattern.
  If \eqn{\gamma > 1}{gamma > 1},
  the model is ``ordered'' or ``inhibitive'' up to the distance
  \eqn{h}, but has an ``attraction'' between points lying at
  distances in the range between \eqn{h} and \eqn{r}.
  
  The function \code{\link{ppm}()}, which fits point process models to 
  point pattern data, requires an argument 
  of class \code{"interact"} describing the interpoint interaction
  structure of the model to be fitted. 
  The appropriate description of the Strauss/hard core process
  pairwise interaction is
  yielded by the function \code{StraussHard()}. See the examples below.
 
  The canonical parameter \eqn{\log(\gamma)}{log(gamma)}
  is estimated by \code{\link{ppm}()}, not fixed in
  \code{StraussHard()}. 
}
\seealso{
  \code{\link{ppm}},
  \code{\link{pairwise.family}},
  \code{\link{ppm.object}}
}
\references{
  Baddeley, A. and Turner, R. (2000)
  Practical maximum pseudolikelihood for spatial point patterns.
  \emph{Australian and New Zealand Journal of Statistics}
  \bold{42}, 283--322.

  Ripley, B.D. (1981) 
  \emph{Spatial statistics}.
  John Wiley and Sons.

  Strauss, D.J. (1975)
  A model for clustering.
  \emph{Biometrika} \bold{63}, 467--475.
}
\examples{
   StraussHard(r=1,hc=0.02)
   # prints a sensible description of itself

   data(cells)

   \dontrun{
   ppm(cells, ~1, StraussHard(r=0.1, hc=0.05))
   # fit the stationary Strauss/hard core  process to `cells'
   }

   ppm(cells, ~ polynom(x,y,3), StraussHard(r=0.1, hc=0.05))
   # fit a nonstationary Strauss/hard core process
   # with log-cubic polynomial trend
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{models}
