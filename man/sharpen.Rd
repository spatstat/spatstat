\name{sharpen} %DontDeclareMethods
\alias{sharpen}
\alias{sharpen.ppp}
\title{Data Sharpening of Point Pattern}
\description{
  Performs Choi-Hall data sharpening of a spatial point pattern.
}
\usage{
sharpen(X, ...)
\method{sharpen}{ppp}(X, sigma=NULL, ..., varcov=NULL,
                        edgecorrect=FALSE)
}
\arguments{
  \item{X}{A marked point pattern (object of class \code{"ppp"}).}
  \item{sigma}{
    Standard deviation of isotropic Gaussian smoothing kernel.
  }
  \item{varcov}{
    Variance-covariance matrix of anisotropic Gaussian kernel.
    Incompatible with \code{sigma}.
  }
  \item{edgecorrect}{
    Logical value indicating whether to apply
    edge effect bias correction.
  }
  \item{\dots}{Arguments passed to \code{\link{density.ppp}}
    to control the pixel resolution of the result.}
}
\details{
  Choi and Hall (2001) proposed a procedure for
  \emph{data sharpening} of spatial point patterns.
  This procedure is appropriate for earthquake epicentres
  and other point patterns which are believed to exhibit
  strong concentrations of points along a curve. Data sharpening
  causes such points to concentrate more tightly along the curve.
  
  If the original data points are 
  \eqn{X_1, \ldots, X_n}{X[1],..., X[n]}
  then the sharpened points are
  \deqn{
    \hat X_i = \frac{\sum_j X_j k(X_j-X_i)}{\sum_j k(X_j - X_i)}
  }{
    X^[i] = (sum[j] X[j] * k(X[j] - X[i]))/(sum[j] k(X[j] - X[i]))
  }
  where \eqn{k} is a smoothing kernel in two dimensions.
  Thus, the new point \eqn{\hat X_i}{X^[i]} is a
  vector average of the nearby points \eqn{X[j]}.

  The function \code{sharpen} is generic. It currently has only one
  method, for two-dimensional point patterns (objects of class
  \code{"ppp"}).

  If \code{sigma} is given, the smoothing kernel is the
  isotropic two-dimensional Gaussian density with standard deviation
  \code{sigma} in each axis. If \code{varcov} is given, the smoothing
  kernel is the Gaussian density with variance-covariance matrix
  \code{varcov}.
  
  The data sharpening procedure tends to cause the point pattern
  to contract away from the boundary of the window. That is,
  points \code{X_i}{X[i]} that lie `quite close to the edge of the window
  of the point pattern tend to be displaced inward. 
  If \code{edgecorrect=TRUE} then the algorithm is modified to
  correct this vector bias. 
}
\value{
  A point pattern (object of class \code{"ppp"}) in the same window
  as the original pattern \code{X}, and with the same marks as \code{X}.
}
\seealso{
  \code{\link{density.ppp}},
  \code{\link{Smooth.ppp}}.
}
\examples{
   data(shapley)
   X <- unmark(shapley)
   \dontshow{
   if(!(interactive())) X <- rthin(X, 0.05)
   }
   Y <- sharpen(X, sigma=0.5)
}
\references{
  Choi, E. and Hall, P. (2001)
  Nonparametric analysis of earthquake point-process data.
  In M. de Gunst, C. Klaassen and A. van der Vaart (eds.)
  \emph{State of the art in probability and statistics:
    Festschrift for Willem R. van Zwet},
  Institute of Mathematical Statistics, Beachwood, Ohio.
  Pages 324--344.
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{nonparametric}
