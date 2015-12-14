\name{rPoissonCluster}
\alias{rPoissonCluster}
\title{Simulate Poisson Cluster Process}
\description{
  Generate a random point pattern, a realisation of the
  general Poisson cluster process.
}
\usage{
 rPoissonCluster(kappa, expand, rcluster, win = owin(c(0,1),c(0,1)),
                 \dots, lmax=NULL, nsim=1, drop=TRUE)
}
\arguments{
  \item{kappa}{
    Intensity of the Poisson process of cluster centres.
    A single positive number, a function, or a pixel image.
  }
  \item{expand}{
    Size of the expansion of the simulation window for generating parent
    points. A single non-negative number.
  }
  \item{rcluster}{
    A function which generates random clusters.
  }
  \item{win}{
    Window in which to simulate the pattern.
    An object of class \code{"owin"}
    or something acceptable to \code{\link{as.owin}}.
  }
  \item{\dots}{
    Arguments passed to \code{rcluster}
  }
  \item{lmax}{
    Optional. Upper bound on the values of \code{kappa}
    when \code{kappa} is a function or pixel image.
  }
  \item{nsim}{Number of simulated realisations to be generated.}
  \item{drop}{
    Logical. If \code{nsim=1} and \code{drop=TRUE} (the default), the
    result will be a point pattern, rather than a list 
    containing a point pattern.
  }
}
\value{
  A point pattern (an object of class \code{"ppp"})
  if \code{nsim=1}, or a list of point patterns if \code{nsim > 1}.
  
  Additionally,  some intermediate results of the simulation are
  returned as attributes of the point pattern: see Details.
}
\details{
  This algorithm generates a realisation of the
  general Poisson cluster process, with the cluster mechanism
  given by the function \code{rcluster}. 

  First, the algorithm
  generates a Poisson point process of ``parent'' points 
  with intensity \code{kappa} in an expanded
  window as explained below.. Here \code{kappa} may be a single
  positive number, a function \code{kappa(x, y)}, or a pixel image
  object of class \code{"im"} (see \code{\link{im.object}}).
  See \code{\link{rpoispp}} for details.
  
  Second, each parent point is replaced by a random cluster of points,
  created by calling the function \code{rcluster}.  These clusters are
  combined together to yield a single point pattern, and the restriction
  of this pattern to the window \code{win} is then returned as the
  result of \code{rPoissonCluster}.

  The expanded window consists of \code{\link{as.rectangle}(win)}
  extended by the amount \code{expand} in each direction. The size of
  the expansion is saved in the attribute \code{"expand"} and may be
  extracted by \code{attr(X, "expand")} where \code{X} is the generated
  point pattern.

  The function \code{rcluster} should expect to be called as
  \code{rcluster(xp[i],yp[i],\dots)} for each parent point at a location
  \code{(xp[i],yp[i])}. The return value of \code{rcluster}
  should be a list with elements
  \code{x,y} which are vectors of equal length giving the absolute
  \eqn{x} and \code{y} coordinates of the points in the cluster.

  If the return value of \code{rcluster} is a point pattern (object of
  class \code{"ppp"}) then it may have marks. The result of
  \code{rPoissonCluster} will then be a marked point pattern.

  If required, the intermediate stages of the simulation (the parents
  and the individual clusters) can also be extracted from
  the return value of \code{rPoissonCluster}
  through the attributes \code{"parents"} and \code{"parentid"}.
  The attribute \code{"parents"} is the point pattern of parent points.
  The attribute \code{"parentid"} is an integer vector specifying
  the parent for each of the points in the simulated pattern.
}
\seealso{
  \code{\link{rpoispp}},
  \code{\link{rMatClust}},
  \code{\link{rThomas}},
  \code{\link{rCauchy}},
  \code{\link{rVarGamma}},
  \code{\link{rNeymanScott}},
  \code{\link{rGaussPoisson}}.
}
\examples{
  # each cluster consist of 10 points in a disc of radius 0.2
  nclust <- function(x0, y0, radius, n) {
              return(runifdisc(n, radius, centre=c(x0, y0)))
            }
  plot(rPoissonCluster(10, 0.2, nclust, radius=0.2, n=5))

  # multitype Neyman-Scott process (each cluster is a multitype process)
  nclust2 <- function(x0, y0, radius, n, types=c("a", "b")) {
     X <- runifdisc(n, radius, centre=c(x0, y0))
     M <- sample(types, n, replace=TRUE)
     marks(X) <- M
     return(X)
  }
  plot(rPoissonCluster(15,0.1,nclust2, radius=0.1, n=5))
}
\author{Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  and Rolf Turner \email{r.turner@auckland.ac.nz}
  
}
\keyword{spatial}
\keyword{datagen}

