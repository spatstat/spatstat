\name{varblock}
\alias{varblock}
\title{
  Estimate Variance of Summary Statistic by Subdivision
}
\description{
  This command estimates the variance of
  any summary statistic (such as the \eqn{K}-function)
  by spatial subdivision of a single point pattern dataset.
}
\usage{
varblock(X, fun = Kest, blocks = quadrats(X, nx = nx, ny = ny), ...,
         nx = 3, ny = nx)
}
\arguments{
  \item{X}{
    Point pattern dataset (object of class \code{"ppp"}).
  }
  \item{fun}{
    Function that computes the summary statistic.
  }
  \item{blocks}{
    Optional. A tessellation that specifies the division of
    the space into blocks.
  }
  \item{\dots}{
    Arguments passed to \code{fun}.
  }
  \item{nx,ny}{
    Optional. Number of rectangular blocks
    in the \eqn{x} and \eqn{y} directions.
    Incompatible with \code{blocks}.
  }
}
\details{
  This command computes an estimate of the variance of
  the summary statistic \code{fun(X)} from a single point pattern
  dataset \code{X} using a subdivision method.
  It can be used to plot \bold{confidence intervals}
  for the true value of a summary function such as the \eqn{K}-function.
  
  The window containing \code{X} is divided into pieces by
  an \code{nx * ny} array of rectangles
  (or is divided into pieces of more general shape,
  according to the argument \code{blocks} if it is present).
  The summary statistic \code{fun} is applied to each of the
  corresponding sub-patterns of \code{X} as described below.
  Then the pointwise
  sample mean, sample variance and sample standard deviation
  of these summary statistics are computed. The two-standard-deviation
  confidence intervals are computed.

  The variance is estimated by equation (4.21) of Diggle (2003, page 52).
  This assumes that the point pattern \code{X} is stationary.
  For further details see Diggle (2003, pp 52--53).
  
  The estimate of the summary statistic
  from each block is computed as follows.
  For most functions \code{fun},
  the estimate from block \code{B}
  is computed by finding the subset of \code{X} consisting of
  points that fall inside \code{B},
  and applying \code{fun} to these points, by calling \code{fun(X[B])}.

  However if \code{fun} is the \eqn{K}-function \code{\link{Kest}},
  or any function which has an argument called \code{domain},
  the estimate for each block \code{B} is computed
  by calling \code{fun(X, domain=B)}. In the case of the
  \eqn{K}-function this means that the estimate from block \code{B}
  is computed by counting pairs of
  points in which the \emph{first} point lies in \code{B},
  while the second point may lie anywhere. 
}
\section{Errors}{
  If the blocks are too small, there may be insufficient data
  in some blocks, and the function \code{fun} may report an error.
  If this happens, you need to take larger blocks.
  
  An error message about incompatibility may occur.
  The different function estimates may be incompatible in some cases,
  for example, because they use different default edge corrections
  (typically because the tiles of the tessellation are not the same kind
  of geometric object as the window of \code{X}, or because the default
  edge correction depends on the number of points). To prevent
  this, specify the choice of edge correction,
  in the \code{correction} argument to \code{fun}, if it has one.
  Some edge corrections are only available if you have set
  \code{spatstat.options(gpclib=TRUE)}.

  An alternative to \code{varblock} is Loh's mark bootstrap
  \code{\link{lohboot}}.
}
\value{
  A function value table (object of class \code{"fv"})
  that contains the result of \code{fun(X)} as well as
  the sample mean, sample variance and sample standard deviation
  of the block estimates, together with 
  the upper and lower two-standard-deviation confidence limits.
}
\references{
  Diggle, P.J. (2003)
  \emph{Statistical analysis of spatial point patterns},
  Second edition. Arnold.
}
\author{
  Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\seealso{
 \code{\link{tess}},
 \code{\link{quadrats}} for basic manipulation.
 
 \code{\link{lohboot}} for an alternative bootstrap technique.
}
\examples{
   v <- varblock(amacrine, Kest, nx=4, ny=2)
   v <- varblock(amacrine, Kcross, nx=4, ny=2)
   if(interactive()) plot(v, iso ~ r, shade=c("hiiso", "loiso"))
}
\keyword{nonparametric}
\keyword{spatial}

