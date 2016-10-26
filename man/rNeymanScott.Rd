\name{rNeymanScott}
\alias{rNeymanScott}
\title{Simulate Neyman-Scott Process}
\description{
  Generate a random point pattern, a realisation of the
  Neyman-Scott cluster process.
}
\usage{
 rNeymanScott(kappa, expand, rcluster, win = owin(c(0,1),c(0,1)),
              \dots, lmax=NULL, nsim=1, drop=TRUE,
              nonempty=TRUE, saveparents=TRUE)
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
    A function which generates random clusters,
    or other data specifying the random cluster mechanism.
    See Details.
  }
  \item{win}{
    Window in which to simulate the pattern.
    An object of class \code{"owin"}
    or something acceptable to \code{\link{as.owin}}.
  }
  \item{\dots}{
    Arguments passed to \code{rcluster}.
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
  \item{nonempty}{
    Logical. If \code{TRUE} (the default), a more efficient algorithm is
    used, in which parents are generated conditionally on having at
    least one offspring point. If \code{FALSE}, parents are generated
    even if they have no offspring. Both choices are valid; the default
    is recommended unless you need to simulate all the parent points
    for some other purpose.
  }
  \item{saveparents}{
    Logical value indicating whether to save the locations of the
    parent points as an attribute.
  }
}
\value{
  A point pattern (an object of class \code{"ppp"}) if \code{nsim=1},
  or a list of point patterns if \code{nsim > 1}.
  
  Additionally,  some intermediate results of the simulation are
  returned as attributes of this point pattern: see Details.
}
\details{
  This algorithm generates a realisation of the
  general Neyman-Scott process, with the cluster mechanism
  given by the function \code{rcluster}. 

  First, the algorithm generates a Poisson point process of
  \dQuote{parent} points with intensity \code{kappa} in an expanded
  window as explained below. Here \code{kappa}
  may be a single positive number,
  a function \code{kappa(x,y)},
  or a pixel image object of class \code{"im"} (see
  \code{\link{im.object}}).  See \code{\link{rpoispp}} for details.

  Second, each parent point is replaced by a random cluster
  of points. These clusters are combined together to yield a
  single point pattern, and the restriction of this pattern to the
  window \code{win} is then returned as the result of
  \code{rNeymanScott}.

  The expanded window consists of \code{\link{as.rectangle}(win)}
  extended by the amount \code{expand} in each direction. The size of
  the expansion is saved in the attribute \code{"expand"} and may be
  extracted by \code{attr(X, "expand")} where \code{X} is the generated
  point pattern.  

  The argument \code{rcluster} specifies the cluster mechanism.
  It may be either:
  \itemize{
    \item
    A \code{function} which will be called to generate each random
    cluster (the offspring points of each parent point).
    The function should expect to be called
    in the form \code{rcluster(x0,y0,\dots)} for a parent point at a location
    \code{(x0,y0)}. The return value of \code{rcluster}
    should specify the coordinates of the points in the cluster;
    it may be a list containing elements
    \code{x,y}, or a point pattern (object of
    class \code{"ppp"}). If it is a marked point pattern then the result of
    \code{rNeymanScott} will be a marked point pattern.
    \item
    A \code{list(mu, f)} where \code{mu} specifies the mean
    number of offspring points in each cluster, and \code{f}
    generates the random displacements (vectors pointing from the parent
    to the offspring). In this case, the number of offspring
    in a cluster is assumed to have a Poisson distribution, implying
    that the Neyman-Scott process is also a Cox process.
    The first element \code{mu} should be either a single nonnegative
    number (interpreted as the mean of the Poisson distribution of
    cluster size)
    or a pixel image or a \code{function(x,y)} giving a spatially
    varying mean cluster size (interpreted in the sense of
    Waagepetersen, 2007).
    The second element \code{f} should be a function that will be
    called once in the form \code{f(n)} to generate \code{n} independent
    and identically distributed displacement vectors (i.e. as if there
    were a cluster of size \code{n} with a parent at the origin
    \code{(0,0)}). 
    The function should return
    a point pattern (object of class \code{"ppp"})
    or something acceptable to \code{\link[grDevices]{xy.coords}}
    that specifies the coordinates of \code{n} points. 
  }

  If required, the intermediate stages of the simulation (the
  parents and the individual clusters) can also be extracted from
  the return value of \code{rNeymanScott} through the attributes
  \code{"parents"} and \code{"parentid"}.  The attribute
  \code{"parents"} is the point pattern of parent points.
  The attribute \code{"parentid"} is an integer vector specifying
  the parent for each of the points in the simulated pattern.

  Neyman-Scott models where \code{kappa} is a single number
  and \code{rcluster = list(mu,f)} can be fitted to data
  using the function \code{\link{kppm}}.
}
\section{Inhomogeneous Neyman-Scott Processes}{
  There are several different ways of specifying a spatially inhomogeneous
  Neyman-Scott process:
  \itemize{
    \item
    The point process of parent points can be inhomogeneous. 
    If the argument \code{kappa} is a \code{function(x,y)} or a pixel
    image (object of class \code{"im"}), then it is taken as specifying
    the intensity function of an inhomogeneous Poisson process according
    to which the parent points are generated.
    \item
    The number of points in a typical cluster can
    be spatially varying.
    If the argument \code{rcluster} is a list of two elements
    \code{mu, f} and the first entry \code{mu} is a 
    \code{function(x,y)} or a pixel image (object of class \code{"im"}),
    then \code{mu} is interpreted as the reference intensity
    for offspring points, in the sense of Waagepetersen (2007).
    For a given parent point, the offspring constitute a Poisson process
    with intensity function equal to \code{mu(x, y) * g(x-x0, y-y0)}
    where \code{g} is the probability density of the offspring
    displacements generated by the function \code{f}.

    Equivalently, clusters are first generated with a constant
    expected number of points per cluster: the constant is \code{mumax}, the
    maximum of \code{mu}. Then the offspring are randomly \emph{thinned}
    (see \code{\link{rthin}}) with spatially-varying retention
    probabilities given by \code{mu/mumax}.  
    \item
    The entire mechanism for generating a cluster can
    be dependent on the location of the parent point.
    If the argument \code{rcluster} is a function,
    then the cluster associated with a parent point at location
    \code{(x0,y0)} will be generated by calling
    \code{rcluster(x0, y0, \dots)}. The behaviour of this function
    could depend on the location \code{(x0,y0)} in any fashion.
  }

  Note that if \code{kappa} is an
  image, the spatial domain covered by this image must be large
  enough to include the \emph{expanded} window in which the parent
  points are to be generated. This requirement means that \code{win} must
  be small enough so that the expansion of \code{as.rectangle(win)}
  is contained in the spatial domain of \code{kappa}.  As a result,
  one may wind up having to simulate the process in a window smaller
  than what is really desired.

  In the first two cases, the intensity of the Neyman-Scott process
  is equal to \code{kappa * mu} if at least one of \code{kappa} or
  \code{mu} is a single number, and is otherwise equal to an
  integral involving \code{kappa}, \code{mu} and \code{f}.
}
\seealso{
  \code{\link{rpoispp}},
  \code{\link{rThomas}},
  \code{\link{rGaussPoisson}},
  \code{\link{rMatClust}},
  \code{\link{rCauchy}},
  \code{\link{rVarGamma}}
}
\examples{
  # each cluster consist of 10 points in a disc of radius 0.2
  nclust <- function(x0, y0, radius, n) {
              return(runifdisc(n, radius, centre=c(x0, y0)))
            }
  plot(rNeymanScott(10, 0.2, nclust, radius=0.2, n=5))

  # multitype Neyman-Scott process (each cluster is a multitype process)
  nclust2 <- function(x0, y0, radius, n, types=c("a", "b")) {
     X <- runifdisc(n, radius, centre=c(x0, y0))
     M <- sample(types, n, replace=TRUE)
     marks(X) <- M
     return(X)
  }
  plot(rNeymanScott(15,0.1,nclust2, radius=0.1, n=5))
}
\references{
  Neyman, J. and Scott, E.L. (1958)
  A statistical approach to problems of cosmology.
  \emph{Journal of the Royal Statistical Society, Series B}
  \bold{20}, 1--43.

  Waagepetersen, R. (2007)
  An estimating function approach to inference for
  inhomogeneous Neyman-Scott processes.
  \emph{Biometrics} \bold{63}, 252--258.
}
\author{\adrian
  
  
  and \rolf
  
}
\keyword{spatial}
\keyword{datagen}

